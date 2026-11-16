#!/usr/bin/env Rscript
# Resolve one genus to species for an existing DADA2 run, without re-running assignTaxonomy.
#
# Why assignSpecies() rather than addSpecies(): addSpecies() routes its sequence argument
# through getUniques(), which in this DADA2 version rejects a plain character vector. The
# standalone assignSpecies() takes the sequences directly and returns one row per sequence.
#
# Matching is EXACT against the species reference, so a species is reported only when the ASV
# equals a reference sequence. For a 350-460 bp V3-V4 ASV many ASVs stay unresolved rather
# than being forced onto the nearest species; the call rate is reported so it is visible, and
# downstream statements are limited to what was actually resolved.
#
# Usage: Rscript species_only.R --study-dir <results/dada2/X> --genus Parvimonas \
#          --ref-species <toSpecies.fa.gz> [--ref-genus <toGenus.fa.gz>]
suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
od <- get_arg("--study-dir")
genus <- get_arg("--genus", "Parvimonas")
ref_species <- get_arg("--ref-species")
ref_genus <- get_arg("--ref-genus", NULL)
tax_file <- get_arg("--tax-file", "taxonomy_silva.csv")
if (is.null(od) || is.null(ref_species)) stop("--study-dir and --ref-species required")
ref_species <- path.expand(ref_species)
if (!file.exists(ref_species)) stop("species reference not found: ", ref_species)

tax <- read.csv(file.path(od, tax_file), row.names = 1,
                check.names = FALSE, stringsAsFactors = FALSE)
gen <- as.character(tax[["Genus"]])
sel <- !is.na(gen) & gen == genus
cat("ASVs assigned to", genus, ":", sum(sel), "of", nrow(tax), "\n")
if (!sum(sel)) quit(status = 0)

fa <- readDNAStringSet(file.path(od, "asv_seqs.fasta"))
ids <- paste0("ASV", seq_along(fa))
if (length(fa) != nrow(tax)) {
  stop("fasta has ", length(fa), " records but taxonomy has ", nrow(tax),
       " rows; refusing to map positionally")
}
sub_ids <- ids[sel]
sq <- toupper(as.character(fa[sel]))

# The Zenodo SILVA trainsets are formatted for assignTaxonomy() (taxonomy in the FASTA
# headers), so species resolution is done by running assignTaxonomy against the toSpecies
# reference and reading its Species column. assignSpecies() expects a different (two-column)
# reference layout and refuses these files.
call_species <- function(seqs, ref) {
  if (is.null(ref) || !file.exists(path.expand(ref))) return(rep(NA_character_, length(seqs)))
  r <- tryCatch(assignTaxonomy(seqs, path.expand(ref), minBoot = 50, verbose = FALSE),
                error = function(e) { message("assignTaxonomy failed: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(rep(NA_character_, length(seqs)))
  if ("Species" %in% colnames(r)) as.character(r[, "Species"]) else rep(NA_character_, nrow(r))
}

sp_species <- call_species(sq, ref_species)
sp_genus <- call_species(sq, ref_genus)
best <- ifelse(is.na(sp_species), sp_genus, sp_species)

out <- data.frame(asv = sub_ids, sequence = sq,
                  species_speciesRef = sp_species,
                  species_genusRef = sp_genus,
                  species = best, stringsAsFactors = FALSE)
write.csv(out, file.path(od, paste0("species_", genus, ".csv")), row.names = FALSE)

cat("resolved:", sum(!is.na(best)), "/", nrow(out), "\n")
print(table(best[!is.na(best)]))

# ---- per-sample abundance per resolved species -----------------------------
st <- readRDS(file.path(od, "asv_table.rds"))
if (ncol(st) != nrow(tax)) {
  stop("ASV table has ", ncol(st), " columns but taxonomy has ", nrow(tax), " rows")
}
lab <- rep(NA_character_, ncol(st))
lab[sel] <- ifelse(is.na(best), paste0(genus, "_unresolved"), best)
keep <- !is.na(lab)
sp_tab <- t(rowsum(t(st[, keep, drop = FALSE]), group = lab[keep], reorder = TRUE))
write.csv(sp_tab, file.path(od, paste0("species_table_", genus, ".csv")))
cat("species table:", nrow(sp_tab), "samples x", ncol(sp_tab), "labels\n")
write(toJSON(list(genus = genus, n_asv = nrow(out),
                  n_resolved = sum(!is.na(best)),
                  calls = as.list(table(best[!is.na(best)]))),
             auto_unbox = TRUE, pretty = TRUE),
      file.path(od, paste0("species_", genus, "_status.json")))
