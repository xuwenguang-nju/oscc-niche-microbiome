#!/usr/bin/env Rscript
# Per-cohort DADA2 processing for the OSCC intratumoral-microbiome study.
#
# One cohort per invocation. Cohorts are NEVER co-processed: the six libraries
# differ in platform (MiSeq / HiSeq / Ion Torrent), read length (251-460 bp) and
# depth (23k-968k), so ASVs are inferred per cohort and compared later at genus
# level through per-cohort effect sizes (see the pipeline design, S1/S5).
#
# Usage
#   Rscript run_dada2.R --study PRJNA1415359
#   Rscript run_dada2.R --study PRJNA866676 --max-files 6
#   Rscript run_dada2.R --all
#
# Outputs per cohort (under results/dada2/<study>/)
#   asv_table.rds           sequence table after chimera removal
#   asv_table.csv           same, ASV ids instead of sequences
#   asv_seqs.fasta          ASV sequences
#   taxonomy.rds/.csv       only when a SILVA training set is present
#   qc.csv                  reads in -> filtered -> denoised -> merged -> nonchim
#   run_status.json         machine-readable success/failure record

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# project root = 4 levels above this script: scripts -> 02-dada2 -> modules -> analysis
script_dir <- tryCatch({
  a <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("--file=", "", a[1]))) else NA
}, error = function(e) NA)
root <- if (!is.na(script_dir) && nzchar(script_dir)) {
  normalizePath(file.path(script_dir, "..", "..", "..", ".."), winslash = "/")
} else {
  normalizePath(getwd(), winslash = "/")
}
if (!dir.exists(file.path(root, "results"))) {
  stop("cannot locate project root; run from the project directory or pass --params")
}

# DADA2 lives in the project-local R library (analysis/r-lib), not the system library
local_lib <- file.path(root, "analysis", "r-lib")
.libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
has_flag <- function(flag) flag %in% args


params_file <- get_arg("--params", file.path(root, "analysis", "modules", "02-dada2",
                                             "input", "cohort_params.yaml"))
raw_dir    <- get_arg("--raw",    file.path(root, "data", "raw", "16s"))
out_root   <- get_arg("--out",    file.path(root, "results", "dada2"))
tax_dir    <- get_arg("--taxonomy", file.path(root, "data", "reference", "silva"))
max_files  <- as.integer(get_arg("--max-files", "0"))
# maxEE override exists so the Ion Torrent cohort can be calibrated empirically
# instead of inheriting the Illumina default of 2
maxee_override <- get_arg("--maxee", NULL)
# trimLeft override exists to settle cohorts with ambiguous primer residue
# (PRJNA666746: primer at 5' in only 34%/47% of reads) by testing both settings
trimleft_override <- get_arg("--trimleft", NULL)
# --downsample N caps reads per sample after filtering (DADA2 runtime scales
# with unique sequences, so over-sequenced cohorts otherwise take days)
downsample_override <- get_arg("--downsample", NULL)
multithread <- !has_flag("--no-multithread")
# --subset <tsv> limits processing to listed runs. Used to restart a large cohort
# on just the paired samples needed downstream, which keeps a single invocation
# short enough to survive in this environment.
subset_file <- get_arg("--subset", NULL)
subset_runs <- NULL
if (!is.null(subset_file)) {
  if (!file.exists(subset_file)) stop("subset file not found: ", subset_file)
  st <- read.delim(subset_file, stringsAsFactors = FALSE)
  subset_runs <- unique(as.character(st$run))
  message("subset: ", length(subset_runs), " runs from ", basename(subset_file))
}

params <- yaml::read_yaml(params_file)$cohorts
studies <- if (has_flag("--all")) names(params) else {
  s <- get_arg("--study")
  if (is.null(s)) stop("give --study <name> or --all")
  s
}

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- helpers ----
list_inputs <- function(study) {
  d <- file.path(raw_dir, study)
  if (!dir.exists(d)) stop("no reads for ", study, " at ", d)
  f1 <- sort(list.files(d, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
  if (length(f1) == 0) {
    # single-end library
    f1 <- sort(list.files(d, pattern = "\\.fastq\\.gz$", full.names = TRUE))
    f1 <- f1[!grepl("_[12]\\.fastq\\.gz$", f1)]
    return(list(R1 = f1, R2 = NULL))
  }
  f2 <- sub("_1\\.fastq\\.gz$", "_2\\.fastq\\.gz", f1)
  keep <- file.exists(f2)
  list(R1 = f1[keep], R2 = f2[keep])
}

sample_name <- function(f) sub("\\.fastq\\.gz$", "", basename(f))

# Downsample filtered FASTQ to at most `n` reads, keeping R1/R2 in step.
#
# Why this exists: DADA2's runtime scales with UNIQUE SEQUENCES, not read count.
# PRJNA666746 has ~240k unique sequences per sample (vs ~16k for PRJNA1415359),
# and the full run consumed 28.8 CPU-hours without passing sample 16 of 100.
# Genus-level composition needs only relative abundance, and 100k reads is far
# more than that requires, so capping depth removes most of the cost without
# touching the analysis unit. The chosen depth is recorded in run_status.json.
subsample_fastq <- function(f1, f2, o1, o2, n, seed = 42L) {
  set.seed(seed)
  r1 <- readFastq(f1)
  total <- length(r1)
  if (total <= n) {
    writeFastq(r1, o1, compress = TRUE)
    removed <- 0L
    if (!is.null(f2)) {
      r2 <- readFastq(f2)
      if (length(r2) != total) {
        stop("R1/R2 read counts differ for ", basename(f1), ": ",
             total, " vs ", length(r2))
      }
      writeFastq(r2, o2, compress = TRUE)
    }
    return(c(total = total, kept = total, removed = 0L))
  }
  idx <- sort(sample.int(total, n))
  writeFastq(r1[idx], o1, compress = TRUE)
  if (!is.null(f2)) {
    r2 <- readFastq(f2)
    if (length(r2) != total) {
      stop("R1/R2 read counts differ for ", basename(f1), ": ",
           total, " vs ", length(r2))
    }
    writeFastq(r2[idx], o2, compress = TRUE)
  }
  c(total = total, kept = n, removed = total - n)
}

run_one <- function(study) {
  p <- params[[study]]
  message("\n=== ", study, "  (", p$platform, ", ", p$region,
          if (isTRUE(p$single_end)) ", single-end" else ", paired-end", ") ===")
  t0 <- Sys.time()

  inp <- list_inputs(study)
  if (!is.null(subset_runs)) {
    keep <- sub("_[12]\\.fastq\\.gz$", "", basename(inp$R1)) %in% subset_runs
    if (!any(keep)) stop("subset selected no files for ", study)
    inp$R1 <- inp$R1[keep]
    if (!is.null(inp$R2)) inp$R2 <- inp$R2[keep]
    message("  subset kept ", length(inp$R1), " runs")
  }
  if (max_files > 0) {
    inp$R1 <- head(inp$R1, max_files)
    if (!is.null(inp$R2)) inp$R2 <- head(inp$R2, max_files)
  }
  message("  input files: ", length(inp$R1),
          if (!is.null(inp$R2)) paste0(" pairs (", length(inp$R2), " R2)"))
  if (length(inp$R1) == 0) stop("no input files")

  od <- file.path(out_root, study)
  fd <- file.path(od, "filtered")
  dir.create(fd, showWarnings = FALSE, recursive = TRUE)

  # ---- optional downsampling of RAW reads (before filtering) ----
  #
  # Ordering matters. Filtering one 928k-read sample took >50 s, so all 200 files
  # cost >2.7 h before denoising even began. Subsampling first cuts that
  # decoding/writing work by the subsampling factor. The raw depth is set so
  # that, after the cohort's measured ~75% filter retention, each sample lands
  # near the requested analysis depth.
  ds_raw <- p$downsample %||% 0
  ds_n <- if (!is.null(downsample_override)) {
    as.integer(downsample_override)
  } else as.numeric(ds_raw)
  downsample_applied <- NA
  if (!is.na(ds_n) && ds_n > 0) {
    RETENTION_GUESS <- 0.75
    raw_depth <- ceiling(ds_n / RETENTION_GUESS)
    rd <- file.path(od, "raw_ds")
    dir.create(rd, showWarnings = FALSE, recursive = TRUE)
    sn_raw <- sample_name(inp$R1)
    dsR1 <- file.path(rd, paste0(sn_raw, "_1.rawds.fastq.gz"))
    dsR2 <- if (!is.null(inp$R2)) {
      file.path(rd, paste0(sn_raw, "_2.rawds.fastq.gz"))
    } else NULL
    message("  downsampling RAW reads to ", raw_depth,
            " per sample (analysis depth after filtering ~", ds_n, ")")
    t_ds <- Sys.time()
    for (i in seq_along(inp$R1)) {
      subsample_fastq(inp$R1[i], if (is.null(dsR2)) NULL else inp$R2[i],
                      dsR1[i], if (is.null(dsR2)) NULL else dsR2[i], raw_depth)
      if (i %% 25 == 0 || i == length(inp$R1)) {
        message("    raw downsampled ", i, "/", length(inp$R1))
      }
    }
    message("  raw downsampling done in ",
            round(as.numeric(difftime(Sys.time(), t_ds, units = "secs")), 1), "s")
    inp$R1 <- dsR1
    if (!is.null(dsR2)) inp$R2 <- dsR2
    downsample_applied <- ds_n
  }

  sn <- sample_name(inp$R1)
  filtR1 <- file.path(fd, paste0(sn, "_F_filt.fastq.gz"))
  filtR2 <- if (!is.null(inp$R2)) file.path(fd, paste0(sn, "_R_filt.fastq.gz")) else NULL

  # ---- filter ----
  # trimLeft comes from the measured primer residue, not from a template.
  tl <- if (!is.null(trimleft_override)) {
    as.integer(strsplit(trimleft_override, ",")[[1]])
  } else as.integer(p$trimLeft %||% 0)
  if (length(tl) == 1) tl <- if (isTRUE(p$single_end)) tl else rep(tl, 2)
  pl_raw <- p$pool %||% FALSE
  pl <- if (is.character(pl_raw)) {
    if (toupper(pl_raw) == "PSEUDO") "pseudo" else as.logical(pl_raw)
  } else as.logical(pl_raw)
  if (is.na(pl)) pl <- FALSE
  message("  pool = ", pl)
  message("  trimLeft = ", paste(tl, collapse = ", "),
          "  truncLen = ", paste(as.integer(p$truncLen), collapse = ", "))
  if (isTRUE(p$single_end)) {
    trunc <- if (length(p$truncLen) > 0 && p$truncLen[1] > 0) as.integer(p$truncLen[1]) else 0L
    me <- if (!is.null(maxee_override)) as.numeric(maxee_override) else as.numeric(p$maxEE %||% 2)
    message("  maxEE = ", me)
    fr <- filterAndTrim(inp$R1, filtR1, truncLen = trunc, trimLeft = tl[1],
                        maxN = 0, maxEE = me, truncQ = 2,
                        minLen = 150, maxLen = 500, rm.phix = FALSE,
                        compress = TRUE, multithread = multithread, verbose = TRUE)
  } else {
    me <- if (!is.null(maxee_override)) as.numeric(maxee_override) else as.numeric(p$maxEE %||% 2)
    message("  maxEE = ", me)
    fr <- filterAndTrim(inp$R1, filtR1, inp$R2, filtR2,
                        truncLen = as.integer(p$truncLen),
                        trimLeft = tl[1:2],
                        maxN = 0, maxEE = c(me, me), truncQ = 2,
                        minLen = 150, maxLen = 500, rm.phix = FALSE,
                        compress = TRUE, multithread = multithread, verbose = TRUE)
  }
  keep <- fr[, "reads.out"] > 0
  if (!any(keep)) stop("every sample had zero reads after filtering")
  filtR1 <- filtR1[keep]
  if (!is.null(filtR2)) filtR2 <- filtR2[keep]
  message("  samples with reads after filtering: ", sum(keep), "/", nrow(fr))

  # ---- learn errors ----
  errF <- learnErrors(filtR1, nbases = 1e8, multithread = multithread, verbose = TRUE)
  if (!is.null(filtR2)) {
    errR <- learnErrors(filtR2, nbases = 1e8, multithread = multithread, verbose = TRUE)
  }

  # ---- denoise ----
  if (isTRUE(p$single_end)) {
    if (isTRUE(p$ion_torrent)) {
      message("  Ion Torrent parameters: HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32")
      dd <- dada(filtR1, err = errF, pool = pl, multithread = multithread,
                 verbose = TRUE, HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32)
    } else {
      dd <- dada(filtR1, err = errF, pool = pl,
                 multithread = multithread, verbose = TRUE)
    }
    seqtab <- makeSequenceTable(dd)
  } else {
    ddF <- dada(filtR1, err = errF, pool = pl,
                multithread = multithread, verbose = TRUE)
    ddR <- dada(filtR2, err = errR, pool = pl,
                multithread = multithread, verbose = TRUE)
    merged <- mergePairs(ddF, filtR1, ddR, filtR2, verbose = TRUE)
    seqtab <- makeSequenceTable(merged)
  }

  # ---- chimera removal ----
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = multithread, verbose = TRUE)

  # ---- taxonomy (optional) ----
  tax <- NULL
  silva_genus <- file.path(tax_dir, "silva_nr99_v138.2_toGenus_trainset.fa.gz")
  if (file.exists(silva_genus)) {
    message("  assigning taxonomy against ", basename(silva_genus))
    tax <- assignTaxonomy(seqtab_nochim, silva_genus, multithread = multithread,
                          taxLevels = c("Kingdom", "Phylum", "Class", "Order",
                                        "Family", "Genus"), verbose = TRUE)
  } else {
    message("  no SILVA training set at ", silva_genus, " -> taxonomy skipped")
  }

  # ---- write ----
  dir.create(od, showWarnings = FALSE, recursive = TRUE)
  saveRDS(seqtab_nochim, file.path(od, "asv_table.rds"))
  if (!is.null(tax)) saveRDS(tax, file.path(od, "taxonomy.rds"))

  asv_seqs <- colnames(seqtab_nochim)
  asv_ids <- paste0("ASV", seq_along(asv_seqs))
  tb <- seqtab_nochim
  colnames(tb) <- asv_ids
  rownames(tb) <- sample_name(rownames(seqtab_nochim))
  write.csv(tb, file.path(od, "asv_table.csv"), row.names = TRUE)
  writeXStringSet(DNAStringSet(asv_seqs), file.path(od, "asv_seqs.fasta"),
                  width = 200)
  if (!is.null(tax)) {
    rownames(tax) <- asv_ids
    write.csv(tax, file.path(od, "taxonomy.csv"), row.names = TRUE)
  }

  read_counts <- if (isTRUE(p$single_end)) {
    rowSums(seqtab)
  } else {
    vapply(merged, function(x) sum(getUniques(x)), numeric(1))
  }
  qc <- data.frame(
    study = study,
    sample = sub("_F_filt\\.fastq\\.gz$", "", basename(filtR1)),
    reads_in = fr[keep, "reads.in"],
    reads_out_filter = fr[keep, "reads.out"],
    reads_denoised = as.integer(read_counts),
    reads_nonchimera = rowSums(seqtab_nochim),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  qc$pct_nonchimera_of_in <- round(100 * qc$reads_nonchimera / qc$reads_in, 1)
  write.csv(qc, file.path(od, "qc.csv"), row.names = FALSE)

  st <- list(
    study = study, status = "ok",
    platform = p$platform, single_end = isTRUE(p$single_end),
    region = p$region,
    truncLen = as.integer(p$truncLen),
    trimLeft = as.integer(tl),
    pool = if (is.logical(pl)) pl else "pseudo",
    ion_torrent = isTRUE(p$ion_torrent),
    downsample_reads = downsample_applied,
    n_samples_in = nrow(fr), n_samples_kept = sum(keep),
    n_asv = length(asv_seqs),
    n_asv_postchimera = ncol(seqtab_nochim),
    median_reads_nonchimera = if (length(read_counts)) median(rowSums(seqtab_nochim)) else 0,
    taxonomy_assigned = !is.null(tax),
    dada2_version = as.character(packageVersion("dada2")),
    elapsed_sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  )
  writeLines(jsonlite::toJSON(st, auto_unbox = TRUE, pretty = TRUE),
             file.path(od, "run_status.json"))
  message("  DONE ", study, ": ", st$n_asv_postchimera, " ASVs, ",
          st$n_samples_kept, " samples, ", st$elapsed_sec, "s")
  st
}

# ---------------------------------------------------------------- driver -----
results <- list()
for (s in studies) {
  results[[s]] <- tryCatch(run_one(s), error = function(e) {
    msg <- conditionMessage(e)
    message("  FAILED ", s, ": ", msg)
    od <- file.path(out_root, s)
    dir.create(od, showWarnings = FALSE, recursive = TRUE)
    writeLines(jsonlite::toJSON(list(study = s, status = "failed", error = msg),
                                auto_unbox = TRUE, pretty = TRUE),
               file.path(od, "run_status.json"))
    list(study = s, status = "failed", error = msg)
  })
}

ok <- vapply(results, function(x) identical(x$status, "ok"), logical(1))
cat("\n==== summary ====\n")
for (s in names(results)) {
  x <- results[[s]]
  cat(sprintf("%-14s %-8s %s\n", s, x$status,
              if (isTRUE(ok[s])) paste0(x$n_asv_postchimera, " ASVs / ",
                                        x$n_samples_kept, " samples") else x$error))
}
if (!all(ok)) quit(status = 1)
