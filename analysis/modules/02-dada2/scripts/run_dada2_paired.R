#!/usr/bin/env Rscript
# Paired-end DADA2 for a cohort whose reads live in <root>/data/raw/16s/<study>.
#
# Self-contained on purpose: the project's run_dada2.R is bound to the project tree and to
# a per-cohort params YAML, which does not exist on a compute host. Parameters here are
# MEASURED from profile_cohort.py output, not assumed:
#   PRJNA412445  2x251 bp MiSeq, no 5' primer residue, Q>=25 to pos 248 (R1) / 245 (R2)
#   -> truncLen 245,240 leaves ~25 bp overlap on a ~460 bp V3-V4 amplicon
#
# Usage: Rscript run_dada2_paired.R --study PRJNA412445 [--trunclen 245,240]
#        [--maxee 2] [--threads 32]
suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
  library(jsonlite)
})

start_time <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
study <- get_arg("--study")
if (is.null(study)) stop("--study is required")
tl_raw <- get_arg("--trunclen", "245,240")
tl <- as.integer(strsplit(tl_raw, ",")[[1]])
maxee <- as.numeric(get_arg("--maxee", "2"))
nthread <- as.integer(get_arg("--threads", "32"))

root <- normalizePath(get_arg("--root", getwd()), winslash = "/")
raw <- file.path(root, "data", "raw", "16s", study)
od <- file.path(root, "results", "dada2", study)
fd <- file.path(od, "filtered")
dir.create(fd, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(root, "dada2_paired.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", file = log_con); flush(log_con)
  message(msg)
}
say("== ", study, " paired-end ==")
say("truncLen=", paste(tl, collapse = ","), " maxEE=", maxee,
    " threads=", nthread)

r1 <- sort(list.files(raw, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
r2 <- sort(list.files(raw, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
if (length(r1) == 0L) stop("no _1 fastq in ", raw)
if (length(r1) != length(r2)) stop("R1/R2 count mismatch")
sn <- sub("_1\\.fastq\\.gz$", "", basename(r1))
say("runs: ", length(r1))

f1 <- file.path(fd, paste0(sn, "_1_filt.fastq.gz"))
f2 <- file.path(fd, paste0(sn, "_2_filt.fastq.gz"))

say("-- filterAndTrim")
fr <- filterAndTrim(r1, f1, r2, f2,
                    truncLen = tl, maxN = 0, maxEE = maxee, truncQ = 2,
                    rm.phix = TRUE, compress = TRUE, multithread = nthread,
                    verbose = FALSE)
say("  pass: ", sum(fr[, 2] > 0), "/", nrow(fr))

keep <- fr[, 2] > 0
if (sum(keep) < 2L) stop("too few samples passed filtering")
k1 <- f1[keep]; k2 <- f2[keep]; kn <- sn[keep]

say("-- learnErrors")
set.seed(42)
e1 <- learnErrors(k1, multithread = nthread, verbose = FALSE)
e2 <- learnErrors(k2, multithread = nthread, verbose = FALSE)
saveRDS(list(forward = e1, reverse = e2), file.path(od, "errors.rds"))

say("-- dada")
d1 <- dada(k1, err = e1, multithread = nthread, pool = FALSE, verbose = FALSE)
d2 <- dada(k2, err = e2, multithread = nthread, pool = FALSE, verbose = FALSE)
say("  R1 denoised; R2 denoised")
saveRDS(list(forward = d1, reverse = d2), file.path(od, "dada.rds"))

say("-- mergePairs")
mrg <- mergePairs(d1, k1, d2, k2, minOverlap = 12, verbose = FALSE)
merged <- makeSequenceTable(mrg)
say("  merged: ", nrow(merged), " samples x ", ncol(merged), " ASVs")

say("-- removeBimeraDenovo")
seqtab <- removeBimeraDenovo(merged, method = "consensus",
                             multithread = nthread, verbose = FALSE)
say("  nonchimera ASVs: ", ncol(seqtab),
    "  reads: ", sum(seqtab))

rownames(seqtab) <- kn
saveRDS(seqtab, file.path(od, "asv_table.rds"))
asv <- colnames(seqtab)
colnames(seqtab) <- paste0("ASV", seq_along(asv))
write.csv(t(seqtab), file.path(od, "asv_table.csv"))
writeLines(asv, file.path(od, "asv_seqs.fasta"))
write.csv(fr, file.path(od, "qc.csv"))

status <- list(
  study = study,
  pipeline = "paired-end",
  trunclen = tl,
  maxee = maxee,
  n_runs = length(r1),
  n_samples_kept = sum(keep),
  n_asv = ncol(seqtab),
  total_reads_nonchimera = sum(seqtab),
  reads_in = sum(fr[, 1]),
  reads_out = sum(fr[, 2]),
  merge_rate = if (nrow(merged) > 0)
    round(sum(getUniques(merged)) / sum(fr[keep, 2]), 4) else NA,
  elapsed_sec = as.numeric(Sys.time() - start_time, units = "secs")
)
write(toJSON(status, auto_unbox = TRUE, pretty = TRUE),
      file.path(od, "run_status.json"))
say("done: ", status$n_asv, " ASVs, ", status$total_reads_nonchimera, " reads")
close(log_con)
