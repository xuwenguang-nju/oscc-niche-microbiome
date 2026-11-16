#!/usr/bin/env Rscript
# Full taxonomic annotation: genus table, family table, and optional species calls.
#
# Three jobs, because each answers a different question:
#   1. genus  -> the tables all downstream analyses (S5/S7) use
#   2. family -> lets us test whether ASVs that fail genus assignment still carry a
#                tumour-enrichment signal at a coarser rank (the "hidden taxa" question)
#   3. species (opt-in) -> resolves a named genus (e.g. Parvimonas) past genus level.
#      Uses addSpecies(), which does exact matching against the species reference, so it
#      reports a species ONLY when the ASV matches a reference sequence exactly. That is
#      deliberately conservative: mis-calling P. micra vs P. dentalis would be worse than
#      returning NA.
#
# Host reads are dropped before collapsing (Eukaryota/Archaea off-target), same rule used
# in the project's assign_taxonomy.R.
#
# Usage: Rscript assign_taxonomy_full.R --study PRJNA866676 --res-root <root>
#        [--species-for Parvimonas] [--threads 16]
suppressPackageStartupMessages({
  library(dada2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
study <- get_arg("--study")
res_root <- get_arg("--res-root", getwd())
ref_root <- get_arg("--ref-root", NULL)
species_for <- get_arg("--species-for", NULL)
nthread <- as.integer(get_arg("--threads", "16"))
if (is.null(study)) stop("--study is required")

if (is.null(ref_root)) ref_root <- file.path(res_root, "data", "reference")
silva_dir <- file.path(ref_root, "silva")
gen_ref <- file.path(silva_dir, "silva_nr99_v138.2_toGenus_trainset.fa.gz")
sp_ref <- file.path(silva_dir, "silva_nr99_v138.2_toSpecies_trainset.fa.gz")

od <- file.path(res_root, "results", "dada2", study)
seqtab <- readRDS(file.path(od, "asv_table.rds"))
seqs <- colnames(seqtab)
cat("study:", study, "| ASVs:", length(seqs), "| samples:", nrow(seqtab), "\n")

RANK6 <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
RANK7 <- c(RANK6, "Species")
tax <- assignTaxonomy(seqs, gen_ref, multithread = nthread, taxLevels = RANK6)
tax_s <- as.data.frame(apply(tax, 2, function(x) {
  y <- as.character(x); y[is.na(y)] <- "unassigned"; y
}), stringsAsFactors = FALSE)
# Do NOT hardcode the number of ranks: the toGenus reference returns 6 (Kingdom..Genus) and
# a species reference would return 7. Renaming by position crashes on the mismatch.
if (!"Genus" %in% colnames(tax_s)) {
  stop("assignTaxonomy returned no Genus column: ", paste(colnames(tax_s), collapse = ","))
}
write.csv(tax_s, file.path(od, "taxonomy_full.csv"))

# ---- host removal, then collapse to genus and to family ---------------------
prok <- tax_s$Kingdom %in% c("Bacteria", "Archaea")
cat("prokaryotic ASVs:", sum(prok), "/", nrow(tax_s), "\n")
st <- seqtab[, prok, drop = FALSE]
tx <- tax_s[prok, , drop = FALSE]

collapse <- function(st, labels, unassigned_label) {
  lab <- labels; lab[lab == "unassigned"] <- unassigned_label
  agg <- rowsum(t(st), group = lab, reorder = TRUE)
  t(agg)
}
genus_tab <- collapse(st, tx$Genus, "unassigned_genus")
family_tab <- collapse(st, tx$Family, "unassigned_family")
write.csv(genus_tab, file.path(od, "genus_table_silva.csv"))
write.csv(family_tab, file.path(od, "family_table_silva.csv"))

fam_na_gen <- tax_s[prok, "Family"] != "unassigned" & tax_s[prok, "Genus"] == "unassigned"
cat("ASVs with family but no genus:", sum(fam_na_gen), "\n")

status <- list(
  study = study,
  n_asv = ncol(seqtab),
  n_prokaryotic = sum(prok),
  n_genera = ncol(genus_tab),
  n_families = ncol(family_tab),
  genus_assigned_pct = round(100 * mean(tax_s$Genus[prok] != "unassigned"), 2),
  reads_total = sum(st),
  genus_table_reads = sum(genus_tab),
  asv_family_no_genus = sum(fam_na_gen)
)

# ---- optional species-level resolution for one genus -----------------------
if (!is.null(species_for) && nzchar(species_for)) {
  target <- tax_s[prok, "Genus"] == species_for
  cat("\nspecies-level for", species_for, ":", sum(target), "ASVs\n")
  if (sum(target) > 0L) {
    sub_seqs <- seqs[prok][target]
    sp <- addSpecies(gen_ref, sub_seqs, allowMultiple = FALSE, verbose = FALSE)
    if (file.exists(sp_ref)) {
      # try the dedicated species reference too; report whichever resolves more exactly
      sp2 <- tryCatch(addSpecies(sp_ref, sub_seqs, allowMultiple = FALSE, verbose = FALSE),
                      error = function(e) NULL)
    } else sp2 <- NULL
    out <- data.frame(
      asv = paste0("ASV", which(prok)[target]),
      sequence = sub_seqs,
      genus = species_for,
      species_genusRef = as.character(sp[, 1]),
      stringsAsFactors = FALSE
    )
    if (!is.null(sp2)) out$species_speciesRef <- as.character(sp2[, 1])
    stt <- st[, target, drop = FALSE]
    out$total_reads <- as.numeric(colSums(stt))
    write.csv(out, file.path(od, paste0("species_", species_for, ".csv")), row.names = FALSE)
    calls <- out$species_genusRef
    if (!is.null(out$species_speciesRef)) {
      calls <- ifelse(is.na(calls), out$species_speciesRef, calls)
    }
    status$species_target <- species_for
    status$species_asv_n <- sum(target)
    status$species_called_n <- sum(!is.na(calls))
    status$species_calls <- as.list(table(calls[!is.na(calls)]))

    # species-resolved abundance matrix for the target genus
    sp_lab <- calls
    sp_lab[is.na(sp_lab)] <- paste0(species_for, "_unresolved")
    sp_tab <- t(rowsum(t(stt), group = sp_lab, reorder = TRUE))
    write.csv(sp_tab, file.path(od, paste0("species_table_", species_for, ".csv")))
    cat("species calls:\n"); print(status$species_calls)
  }
}

write(toJSON(status, auto_unbox = TRUE, pretty = TRUE),
      file.path(od, "taxonomy_status.json"))
cat("\ndone\n")
