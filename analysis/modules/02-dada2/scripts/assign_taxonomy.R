#!/usr/bin/env Rscript
# Assign taxonomy to an existing ASV table (no DADA2 re-run needed).
#
# This exists because taxonomy assignment is the rate-limiting prerequisite for
# the genus-level analysis: every downstream step (decontamination, source
# attribution, cross-cohort meta-analysis) is genus-level. Re-running the whole
# DADA2 pipeline just to attach taxonomy would waste hours.
#
# Two references are applied in sequence and both are kept, because they answer
# different questions:
#   SILVA 138.2  cross-cohort default, consistent across all six cohorts
#   eHOMD V15.22 oral-specific; a disagreement with SILVA is itself informative
#
# Usage
#   Rscript assign_taxonomy.R --study PRJNA1415359
#   Rscript assign_taxonomy.R --all
#
# Outputs per cohort (added to results/dada2/<study>/)
#   taxonomy_silva.rds / .csv
#   taxonomy_ehomd.rds / .csv
#   taxonomy_compare.csv        per-ASV SILVA vs eHOMD genus agreement
#   genus_table_silva.csv       ASV counts collapsed to genus (the analysis unit)
#   taxonomy_status.json

suppressPackageStartupMessages({
  library(yaml); library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

script_dir <- tryCatch({
  a <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("--file=", "", a[1]))) else NA
}, error = function(e) NA)
root <- if (!is.na(script_dir) && nzchar(script_dir)) {
  normalizePath(file.path(script_dir, "..", "..", "..", ".."), winslash = "/")
} else normalizePath(getwd(), winslash = "/")

.libPaths(c(file.path(root, "analysis", "r-lib"), .libPaths()))
suppressPackageStartupMessages({ library(dada2); library(Biostrings) })

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(f) f %in% args
get_arg <- function(f, d = NULL) {
  i <- match(f, args); if (is.na(i) || i == length(args)) return(d); args[i + 1]
}

# --res-root lets this same script run on a compute host where the results tree
# is not laid out under the project root (root is derived from the script path).
res_root <- get_arg("--res-root", file.path(root, "results", "dada2"))
silva <- get_arg("--silva", file.path(root, "data", "reference", "silva",
                                      "silva_nr99_v138.2_toGenus_trainset.fa.gz"))
ehomd <- get_arg("--ehomd", file.path(root, "data", "reference", "ehomd",
                                      "eHOMD_RefSeq_dada2_V15.22.fasta.gz"))
mt <- !has_flag("--no-multithread")

if (!file.exists(silva)) stop("SILVA training set missing: ", silva)
have_ehomd <- file.exists(ehomd)

studies <- if (has_flag("--all")) {
  d <- list.dirs(res_root, recursive = FALSE, full.names = FALSE)
  d[file.exists(file.path(res_root, d, "asv_table.rds"))]
} else {
  s <- get_arg("--study"); if (is.null(s)) stop("give --study or --all"); s
}

collapse_to_genus <- function(seqtab, tax, level = "Genus") {
  # Drop host / non-prokaryotic reads BEFORE collapsing.
  #
  # Tissue-derived 16S libraries routinely carry eukaryotic (host) amplicons that
  # survive filtering, denoising and chimera removal. Measured host read share
  # across these cohorts: PRJNA1415359 0.0%, PRJNA597251 0.0%, PRJNA1415268 2.4%,
  # PRJNA666746 49.0%. Leaving them in dilutes every bacterial genus by that
  # factor, which would bias cross-cohort effect sizes -- so keep only Bacteria
  # and Archaea.
  kingdom <- tax[, "Kingdom"]
  keep <- !is.na(kingdom) & kingdom %in% c("Bacteria", "Archaea")
  seqtab <- seqtab[, keep, drop = FALSE]
  tax <- tax[keep, , drop = FALSE]
  g <- tax[, level]
  g[is.na(g) | g == ""] <- "unassigned"
  # Bacterial ASVs without a genus call stay as "unassigned", so the bacterial
  # read total is conserved and can be checked against reads_prokaryote.
  agg <- rowsum(t(seqtab), group = g, reorder = TRUE)
  t(agg)
}

for (s in studies) {
  od <- file.path(res_root, s)
  st_file <- file.path(od, "taxonomy_status.json")
  message("\n=== ", s, " ===")
  t0 <- Sys.time()
  status <- list(study = s, status = "ok")

  seqtab <- readRDS(file.path(od, "asv_table.rds"))
  seqs <- colnames(seqtab)
  asv_ids <- paste0("ASV", seq_along(seqs))
  colnames(seqtab) <- asv_ids
  message("  ", ncol(seqtab), " ASVs x ", nrow(seqtab), " samples")

  message("  assigning with SILVA 138.2 ...")
  tax_s <- assignTaxonomy(seqs, silva, multithread = mt,
                          taxLevels = c("Kingdom", "Phylum", "Class", "Order",
                                        "Family", "Genus"), verbose = TRUE)
  rownames(tax_s) <- asv_ids
  saveRDS(tax_s, file.path(od, "taxonomy_silva.rds"))
  write.csv(tax_s, file.path(od, "taxonomy_silva.csv"))
  genus_tab <- collapse_to_genus(seqtab, tax_s)
  write.csv(genus_tab, file.path(od, "genus_table_silva.csv"))

  status$n_asv <- ncol(seqtab)
  status$genus_assigned <- sum(!is.na(tax_s[, "Genus"]))
  status$genus_assigned_pct <- round(100 * status$genus_assigned / ncol(seqtab), 1)
  status$n_genera <- ncol(genus_tab)   # genus_tab is samples x genera
  status$reads_total <- sum(seqtab)
  status$reads_prokaryote <- sum(genus_tab)
  # Was previously reads_genus_assigned_pct, which was always 100 because
  # "unassigned" counted as a genus. Report host removal explicitly instead.
  status$reads_prokaryote_pct <-
    round(100 * status$reads_prokaryote / status$reads_total, 1)
  kingdom <- tax_s[, "Kingdom"]
  status$host_reads_pct <-
    round(100 * sum(seqtab[, is.na(kingdom) | !kingdom %in% c("Bacteria", "Archaea")]) /
            status$reads_total, 1)
  status$bacterial_no_genus_reads_pct <-
    round(100 * sum(genus_tab[, "unassigned"]) / status$reads_total, 1)

  if (have_ehomd) {
    message("  assigning with eHOMD V15.22 ...")
    tax_e <- tryCatch(
      assignTaxonomy(seqs, ehomd, multithread = mt,
                     taxLevels = c("Kingdom", "Phylum", "Class", "Order",
                                   "Family", "Genus"), verbose = TRUE),
      error = function(e) { message("    eHOMD failed: ", conditionMessage(e)); NULL })
    if (!is.null(tax_e)) {
      rownames(tax_e) <- asv_ids
      saveRDS(tax_e, file.path(od, "taxonomy_ehomd.rds"))
      write.csv(tax_e, file.path(od, "taxonomy_ehomd.csv"))
      cmp <- data.frame(
        asv = asv_ids,
        genus_silva = tax_s[, "Genus"],
        genus_ehomd = tax_e[, "Genus"],
        row.names = NULL, stringsAsFactors = FALSE)
      cmp$agree <- !is.na(cmp$genus_silva) & !is.na(cmp$genus_ehomd) &
        cmp$genus_silva == cmp$genus_ehomd
      cmp$silva_only <- !is.na(cmp$genus_silva) & is.na(cmp$genus_ehomd)
      cmp$ehomd_only <- is.na(cmp$genus_silva) & !is.na(cmp$genus_ehomd)
      write.csv(cmp, file.path(od, "taxonomy_compare.csv"), row.names = FALSE)
      status$ehomd_genus_assigned <- sum(!is.na(tax_e[, "Genus"]))
      status$ehomd_genus_assigned_pct <-
        round(100 * status$ehomd_genus_assigned / ncol(seqtab), 1)
      status$genus_agreement_pct <-
        round(100 * sum(cmp$agree) / ncol(seqtab), 1)
      # reads-weighted agreement matters more than ASV-count agreement
      w <- colSums(seqtab)   # seqtab is samples x ASVs, so weight per ASV is a column sum
      status$genus_agreement_reads_pct <-
        round(100 * sum(w[cmp$agree]) / sum(w), 1)
    }
  }

  status$elapsed_sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  writeLines(jsonlite::toJSON(status, auto_unbox = TRUE, pretty = TRUE), st_file)
  message("  genus assigned: ", status$genus_assigned_pct, "% of ASVs, ",
          status$reads_genus_assigned_pct, "% of reads; ",
          status$n_genera, " genera")
  if (!is.null(status$genus_agreement_pct)) {
    message("  SILVA vs eHOMD genus agreement: ", status$genus_agreement_pct,
            "% of ASVs, ", status$genus_agreement_reads_pct, "% of reads")
  }
  message("  done in ", status$elapsed_sec, "s")
}
