#!/usr/bin/env Rscript
# Install DADA2 into a PROJECT-LOCAL R library on the Linux compute host.
#
# Why not the system library: /opt/R/4.4.0/lib/R/library is not writable and we
# are unprivileged (no sudo), so every package goes to $OSCC_R_LIB.
#
# Why Bioconductor 3.20 rather than the release that matches R 4.4.0 (3.19):
# the three cohorts already processed locally used DADA2 1.34.0, and using the
# same version on the fourth cohort keeps the pipeline identical across cohorts.
# Bioc 3.20 requires R >= 4.4.0, which this host satisfies.
#
# Needs: gcc/g++/gfortran/make (DADA2 compiles C++). No cmake required.

lib <- Sys.getenv("OSCC_R_LIB", unset = file.path(Sys.getenv("HOME"), "oscc-mb", "Rlib"))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cat("R:", R.version.string, "\n")
cat("lib:", lib, "\n")
cat("libPaths:", paste(.libPaths(), collapse = " : "), "\n")

repos <- c(CRAN = "https://cloud.r-project.org")
options(repos = repos, Ncpus = max(1L, min(16L, parallel::detectCores() - 2L)))
cat("Ncpus:", getOption("Ncpus"), "\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib, repos = repos)
}
# Pin Bioc 3.20 so that dada2 resolves to 1.34.x (matches the local cohorts).
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
cat("BiocManager:", as.character(packageVersion("BiocManager")), "\n")

# jsonlite: run_dada2.R writes run_status.json with it, and this host lacks it.
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", lib = lib, repos = repos)
}
cat("jsonlite:", as.character(packageVersion("jsonlite")), "\n")

BiocManager::install("dada2", lib = lib, ask = FALSE, update = FALSE)

ok <- requireNamespace("dada2", quietly = TRUE)
cat("\n==== RESULT ====\n")
cat("dada2 installed:", ok, "\n")
if (ok) {
  cat("dada2 version:", as.character(packageVersion("dada2")), "\n")
  cat("library path:", dirname(system.file(package = "dada2")), "\n")
  suppressPackageStartupMessages(library(dada2))
  cat("dada() has BAND_SIZE:", "BAND_SIZE" %in% names(formals(dada2::dada)), "\n")
  cat("cores visible:", parallel::detectCores(), "\n")
}
if (!ok) quit(status = 1)
