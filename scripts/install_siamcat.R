local_library <- Sys.getenv(
  "R_LIBS_USER",
  unset = Sys.getenv("IBD_R_LIB", unset = "")
)
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .libPaths()))
message("Using local R library: ", local_library)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
message("Active Bioconductor version: ", BiocManager::version())
message("Installing SIAMCAT and dependencies")
BiocManager::install("SIAMCAT", ask = FALSE, update = FALSE)

if (!requireNamespace("SIAMCAT", quietly = TRUE)) {
  stop("SIAMCAT installation did not complete")
}
cat("SIAMCAT version: ", as.character(utils::packageVersion("SIAMCAT")), "\n", sep = "")
cat("Library: ", find.package("SIAMCAT"), "\n", sep = "")

