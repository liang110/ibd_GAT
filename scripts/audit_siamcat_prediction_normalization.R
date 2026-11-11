#!/usr/bin/env Rscript
# Print the installed SIAMCAT implementation needed to audit target handling.

if (!requireNamespace("SIAMCAT", quietly = TRUE)) {
  stop("SIAMCAT must be installed before this audit can run.")
}

for (function_name in c(
  "make.predictions", "make.external.predictions", "normalize.features"
)) {
  cat("\n###", function_name, "\n")
  print(getAnywhere(function_name))
}

