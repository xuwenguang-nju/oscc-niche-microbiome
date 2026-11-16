# Install DADA2 into a project-local R library (does not touch the system R library).
#
# Repositories: official Bioconductor for Bioc packages (the Tsinghua Bioconductor
# mirror no longer serves /packages/<ver>/... so it 404s on PACKAGES and silently
# breaks resolution), Tsinghua CRAN for CRAN packages.
#
# R 4.4.2 -> Bioconductor 3.20 -> dada2 1.34.0 (prebuilt Windows binary; no Rtools needed).

lib <- file.path("analysis", "r-lib")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
lib <- normalizePath(lib, winslash = "/")
.libPaths(c(lib, .libPaths()))

options(
  repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
  timeout = 1800
)

cat("R:", R.version.string, "\n")
cat("BiocManager:", as.character(packageVersion("BiocManager")), "\n")
cat("target lib:", lib, "\n")
cat("Bioc version:", as.character(BiocManager::version()), "\n")

if (!requireNamespace("dada2", quietly = TRUE)) {
  cat("installing dada2 from Bioconductor 3.20 (official repo) ...\n")
  BiocManager::install("dada2", lib = lib, update = FALSE, ask = FALSE,
                       site_repository = NULL,
                       Ncpus = max(1, parallel::detectCores() - 1))
}

ok <- requireNamespace("dada2", quietly = TRUE)
cat("dada2 available:", ok, "\n")
if (ok) {
  cat("dada2 version:", as.character(packageVersion("dada2")), "\n")
  for (p in c("Biostrings", "ShortRead", "RcppParallel", "IRanges", "S4Vectors")) {
    v <- tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")
    cat(" ", p, v, "\n")
  }
  cat("library path:", dirname(system.file(package = "dada2")), "\n")
  cat("dada() args include HOMOPOLYMER_GAP_PENALTY/BAND_SIZE:",
      all(c("HOMOPOLYMER_GAP_PENALTY", "BAND_SIZE") %in% names(formals(dada2::dada))), "\n")
}
