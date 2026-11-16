#!/usr/bin/env Rscript
# DADA2 for PRJNA866676 — Ion Torrent, single-end, six ecological niches.
#
# WHY A DEDICATED SCRIPT
#   The shared run_dada2.R assumes a project tree plus a multi-cohort YAML; this cohort
#   runs on the compute host from a flat work directory, so the paths and the
#   single-end/Ion-Torrent choices are made explicit here instead of being inferred.
#
# PARAMETERS ARE MEASURED, NOT ASSUMED (see profile_866676.py output)
#   read length : median 426 bp, IQR 419-429, max 460
#   mean quality: Q>=25 up to position 349; Q>=20 up to 407
#   -> truncLen = 350 keeps 100% of reads and 82.8% of bases at Q>=25
#   -> maxEE = 3   (Ion Torrent needs a looser EE than MiSeq's 2; calibrated earlier on
#                   this cohort: maxEE 2->3 roughly doubled retention with only ~9% more
#                   ASVs, i.e. the extra reads were duplicates of detected sequences)
#   -> trimLeft = 0 (no primer residue detected in these reads)
#   -> pool = FALSE (consistent with the other cohorts)
#
# Usage:
#   R_LIBS_USER=~/oscc-mb/Rlib Rscript run_dada2_single_end.R \
#       --raw  data/raw/16s/PRJNA866676 \
#       --out  results/dada2/PRJNA866676 \
#       --trunclen 350 --maxee 3 --threads 32

suppressPackageStartupMessages({
  library(dada2)
  library(jsonlite)
})

get_arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, a)
  if (is.na(i) || i == length(a)) return(default)
  a[i + 1]
}

raw_dir <- get_arg("--raw")
out_dir <- get_arg("--out")
truncLen <- as.integer(get_arg("--trunclen", "350"))
maxEE <- as.numeric(get_arg("--maxee", "3"))
threads <- as.integer(get_arg("--threads", "32"))
study <- get_arg("--study", basename(out_dir))
if (is.null(raw_dir) || is.null(out_dir)) stop("--raw and --out are required")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_con <- file(file.path(out_dir, "dada2_single_end.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", file = log_con); flush(log_con); message(msg)
}

t0 <- Sys.time()
say("== DADA2 single-end (Ion Torrent) ==")
say("host     : ", Sys.info()[["nodename"]], " (", Sys.info()[["sysname"]], ")")
say("dada2    : ", as.character(packageVersion("dada2")))
say("R        : ", R.version.string)
say("threads  : ", threads, " of ", parallel::detectCores())
say("truncLen : ", truncLen, "   maxEE: ", maxEE, "   pool: FALSE")
say("raw      : ", raw_dir)

fnFs <- sort(list.files(raw_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE))
if (length(fnFs) == 0) stop("no fastq.gz in ", raw_dir)
sample_names <- sub("\\.fastq\\.gz$", "", basename(fnFs))
names(fnFs) <- sample_names
say("samples  : ", length(fnFs))

# ---- filter + trim -------------------------------------------------------
say("\n[1/3] filterAndTrim ...")
filt_dir <- file.path(out_dir, "filtered")
dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)
fnFs_filt <- file.path(filt_dir, paste0(sample_names, ".filt.fastq.gz"))

fr <- filterAndTrim(fnFs, fnFs_filt, truncLen = truncLen, maxEE = maxEE,
                    truncQ = 2, maxN = 0, rm.phix = FALSE,   # Ion Torrent: no PhiX
                    compress = TRUE, multithread = TRUE, verbose = TRUE)
saveRDS(fr, file.path(out_dir, "filter_stats.rds"))
write.csv(fr, file.path(out_dir, "filter_stats.csv"))
say("  median retention: ",
    round(100 * median(fr[, 2] / fr[, 1]), 1), "%")
keep <- fr[, 2] > 0
say("  samples with reads: ", sum(keep), "/", nrow(fr))
if (sum(keep) < 10) stop("too few samples retained")

filt <- fnFs_filt[keep]
names(filt) <- sample_names[keep]

# ---- learn errors + denoise ---------------------------------------------
say("\n[2/3] learnErrors + dada ...")
set.seed(42)
err <- learnErrors(filt, multithread = TRUE, nbases = 2e8, verbose = 2)
saveRDS(err, file.path(out_dir, "errors.rds"))

# Ion Torrent needs BAND_SIZE raised because indels shift the alignment band.
dada_out <- dada(filt, err = err, multithread = TRUE, pool = FALSE,
                 BAND_SIZE = 32, verbose = 1)
saveRDS(dada_out, file.path(out_dir, "dada.rds"))
say("  denoised; elapsed ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

seqtab <- makeSequenceTable(dada_out)
say("  table: ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs")

# ---- chimera + save -----------------------------------------------------
say("\n[3/3] removeBimeraDenovo ...")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
say("  non-chimeric: ", nrow(seqtab_nochim), " x ", ncol(seqtab_nochim))
say("  reads retained: ",
    round(100 * sum(seqtab_nochim) / sum(seqtab), 2), "% of denoised")

saveRDS(seqtab_nochim, file.path(out_dir, "asv_table.rds"))
write.csv(t(seqtab_nochim), file.path(out_dir, "asv_table.csv"))
seqs <- colnames(seqtab_nochim)
names(seqs) <- paste0("ASV", seq_along(seqs))
writeLines(paste0(">", names(seqs), "\n", seqs),
           file.path(out_dir, "asv_seqs.fasta"))

qc <- data.frame(
  sample = rownames(seqtab_nochim),
  reads_in = as.integer(fr[rownames(seqtab_nochim), 1]),
  reads_filt = as.integer(fr[rownames(seqtab_nochim), 2]),
  reads_denoised = as.integer(rowSums(seqtab)),
  reads_nonchimera = as.integer(rowSums(seqtab_nochim)),
  asvs = as.integer(rowSums(seqtab_nochim > 0)),
  pct_nonchimera_of_filtered =
    round(100 * rowSums(seqtab_nochim) / fr[rownames(seqtab_nochim), 2], 2),
  row.names = NULL)
write.csv(qc, file.path(out_dir, "qc.csv"), row.names = FALSE)

status <- list(
  study = study, finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  dada2_version = as.character(packageVersion("dada2")),
  r_version = R.version.string, host = Sys.info()[["nodename"]],
  platform = Sys.info()[["sysname"]], threads = threads,
  single_end = TRUE, ion_torrent = TRUE, pool = FALSE, band_size = 32,
  truncLen = truncLen, maxEE = maxEE, truncQ = 2,
  n_samples_in = nrow(fr), n_samples_kept = sum(keep),
  n_asv = ncol(seqtab_nochim),
  total_reads_nonchimera = as.numeric(sum(seqtab_nochim)),
  median_reads_nonchimera = median(rowSums(seqtab_nochim)),
  median_pct_nonchimera_of_filtered = median(qc$pct_nonchimera_of_filtered),
  elapsed_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
write(toJSON(status, auto_unbox = TRUE, pretty = TRUE),
      file.path(out_dir, "run_status.json"))

say("\n== DONE ==")
for (k in names(status)) say("  ", k, ": ", status[[k]])
close(log_con)
