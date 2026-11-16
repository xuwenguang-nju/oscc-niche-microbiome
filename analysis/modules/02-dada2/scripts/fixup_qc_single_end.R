#!/usr/bin/env Rscript
# Rebuild qc.csv and run_status.json for PRJNA866676 without re-running DADA2.
#
# WHY: the DADA2 run itself completed (35,375 ASVs, 17.29 M reads) and all primary
# objects were written, but the QC table aborted because filterAndTrim's rownames are
# input paths while the ASV table is indexed by sample name. That mismatch is fixed in
# run_dada2_single_end.R; this script only regenerates the two derived files.

out_dir <- "results/dada2/PRJNA866676"
suppressPackageStartupMessages(library(jsonlite))

fr <- readRDS(file.path(out_dir, "filter_stats.rds"))
rownames(fr) <- sub("\\.fastq\\.gz$", "", rownames(fr))
seqtab <- readRDS(file.path(out_dir, "asv_table.rds"))
stopifnot(all(rownames(seqtab) %in% rownames(fr)))

qc <- data.frame(
  sample = rownames(seqtab),
  reads_in = as.integer(fr[rownames(seqtab), "reads.in"]),
  reads_filt = as.integer(fr[rownames(seqtab), "reads.out"]),
  reads_nonchimera = as.integer(rowSums(seqtab)),
  asvs = as.integer(rowSums(seqtab > 0)),
  pct_nonchimera_of_filtered =
    round(100 * rowSums(seqtab) / fr[rownames(seqtab), "reads.out"], 2),
  row.names = NULL)
write.csv(qc, file.path(out_dir, "qc.csv"), row.names = FALSE)

status <- list(
  study = "PRJNA866676",
  finished = format(file.mtime(file.path(out_dir, "asv_table.rds")),
                    "%Y-%m-%dT%H:%M:%S%z"),
  dada2_version = as.character(packageVersion("dada2")),
  r_version = R.version.string,
  single_end = TRUE, ion_torrent = TRUE, pool = FALSE, band_size = 32,
  truncLen = 350, maxEE = 3, truncQ = 2,
  n_samples_in = nrow(fr), n_samples_kept = nrow(seqtab),
  n_asv = ncol(seqtab),
  total_reads_nonchimera = as.numeric(sum(seqtab)),
  median_reads_in = median(qc$reads_in),
  median_reads_nonchimera = median(qc$reads_nonchimera),
  median_pct_nonchimera_of_filtered = median(qc$pct_nonchimera_of_filtered),
  note = "qc/status regenerated after the original QC step hit a rowname mismatch; \
the ASV table and sequences are from the original run")
write(toJSON(status, auto_unbox = TRUE, pretty = TRUE),
      file.path(out_dir, "run_status.json"))

cat("qc rows:", nrow(qc), "\n")
cat("median reads in:", status$median_reads_in,
    " median non-chimera:", status$median_reads_nonchimera,
    " median retention:", status$median_pct_nonchimera_of_filtered, "%\n")
cat("ASVs:", status$n_asv, " total reads:", status$total_reads_nonchimera, "\n")
