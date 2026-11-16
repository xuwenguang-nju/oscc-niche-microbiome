#!/usr/bin/env Rscript
# Run DADA2 denoising on an ALREADY-FILTERED set of FASTQ files.
#
# Why this exists instead of reusing run_dada2.R: on this Linux host the
# filterAndTrim step has already been done (on a different machine), so this
# script starts at learnErrors(). Its whole reason to exist is that DADA2's
# multithreading needs fork(), which Windows lacks -- the same step that ran at
# ~1-2 cores locally uses dozens of cores here.
#
# Input naming (written by run_dada2.R after --downsample):
#   <sample>.rawds_F_filt.fastq.gz   forward
#   <sample>.rawds_R_filt.fastq.gz   reverse
#
# Usage:
#   Rscript run_dada2_from_filtered.R --filtered DIR --out DIR [--threads 32]

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

filtered_dir <- get_arg("--filtered")
out_dir <- get_arg("--out")
threads <- as.integer(get_arg("--threads", "32"))
study <- get_arg("--study", basename(normalizePath(out_dir)))
if (is.null(filtered_dir) || is.null(out_dir)) {
  stop("--filtered and --out are required")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_con <- file(file.path(out_dir, "dada2_from_filtered.log"), open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", file = log_con)
  flush(log_con)
  message(msg)
}

t0 <- Sys.time()
say("== DADA2 from filtered FASTQ ==")
say("host            : ", Sys.info()[["nodename"]], " (", Sys.info()[["sysname"]], ")")
say("dada2           : ", as.character(packageVersion("dada2")))
say("R               : ", R.version.string)
say("threads         : ", threads, " of ", parallel::detectCores(), " visible")
say("filtered dir    : ", filtered_dir)
say("out dir         : ", out_dir)

# ---- discover paired filtered files -------------------------------------
fnFs <- sort(list.files(filtered_dir, pattern = "\\.rawds_F_filt\\.fastq\\.gz$",
                        full.names = TRUE))
fnRs <- sort(list.files(filtered_dir, pattern = "\\.rawds_R_filt\\.fastq\\.gz$",
                        full.names = TRUE))
say("forward files   : ", length(fnFs))
say("reverse files   : ", length(fnRs))
if (length(fnFs) == 0) stop("no *_F_filt.fastq.gz found in ", filtered_dir)

sample_names <- sub("\\.rawds_F_filt\\.fastq\\.gz$", "", basename(fnFs))
names(fnFs) <- sample_names
names(fnRs) <- sub("\\.rawds_R_filt\\.fastq\\.gz$", "", basename(fnRs))
stopifnot(identical(sort(names(fnFs)), sort(names(fnRs))))
fnRs <- fnRs[names(fnFs)]

say("paired samples  : ", length(fnFs))

# ---- error learning -----------------------------------------------------
say("\n[1/4] learnErrors ...")
errF <- learnErrors(fnFs, multithread = TRUE, nbases = 2e8, verbose = 2)
errR <- learnErrors(fnRs, multithread = TRUE, nbases = 2e8, verbose = 2)
saveRDS(list(F = errF, R = errR), file.path(out_dir, "errors.rds"))
say("  error model saved; elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

# ---- denoise ------------------------------------------------------------
# pool = FALSE: matches the other three cohorts. These are high-diversity
# tissue/saliva samples, so pseudo-pooling would only inflate runtime.
say("\n[2/4] dada (pool = FALSE, multithread = TRUE) ...")
dadaFs <- dada(fnFs, err = errF, multithread = TRUE, pool = FALSE, verbose = 1)
dadaRs <- dada(fnRs, err = errR, multithread = TRUE, pool = FALSE, verbose = 1)
saveRDS(list(F = dadaFs, R = dadaRs), file.path(out_dir, "dada.rds"))
say("  denoised; elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

# ---- merge --------------------------------------------------------------
say("\n[3/4] mergePairs ...")
mergers <- mergePairs(dadaFs, fnFs, dadaRs, fnRs, verbose = TRUE)
seqtab <- makeSequenceTable(mergers)
say("  merged table    : ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs")

# ---- chimera removal ----------------------------------------------------
say("\n[4/4] removeBimeraDenovo ...")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
say("  non-chimeric    : ", nrow(seqtab_nochim), " samples x ", ncol(seqtab_nochim), " ASVs")
say("  reads retained  : ", round(100 * sum(seqtab_nochim) / sum(seqtab), 2), "% of merged")

saveRDS(seqtab_nochim, file.path(out_dir, "seqtab_nochim.rds"))
write.csv(t(seqtab_nochim), file.path(out_dir, "seqtab_nochim.csv"))

# ---- QC table -----------------------------------------------------------
qc <- data.frame(
  sample = rownames(seqtab_nochim),
  reads_in = as.integer(rowSums(seqtab)),
  reads_nonchimera = as.integer(rowSums(seqtab_nochim)),
  asvs = as.integer(rowSums(seqtab_nochim > 0)),
  pct_nonchimera_of_merged = round(100 * rowSums(seqtab_nochim) / rowSums(seqtab), 2),
  row.names = NULL
)
write.csv(qc, file.path(out_dir, "qc.csv"), row.names = FALSE)

status <- list(
  study = study,
  finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  dada2_version = as.character(packageVersion("dada2")),
  r_version = R.version.string,
  host = Sys.info()[["nodename"]],
  platform = Sys.info()[["sysname"]],
  threads = threads,
  pool = FALSE,
  n_samples = nrow(seqtab_nochim),
  n_asvs = ncol(seqtab_nochim),
  total_reads_nonchimera = as.numeric(sum(seqtab_nochim)),
  median_reads_nonchimera = median(rowSums(seqtab_nochim)),
  median_pct_nonchimera_of_merged = median(qc$pct_nonchimera_of_merged),
  elapsed_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
)
write(toJSON(status, auto_unbox = TRUE, pretty = TRUE),
      file.path(out_dir, "run_status.json"))

say("\n== DONE ==")
for (k in names(status)) say("  ", k, ": ", status[[k]])
close(log_con)
