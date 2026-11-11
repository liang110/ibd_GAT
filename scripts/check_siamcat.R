cat("R version: ", R.version.string, "\n", sep = "")
cat("SIAMCAT installed: ", requireNamespace("SIAMCAT", quietly = TRUE), "\n", sep = "")
cat("R libraries:\n")
writeLines(.libPaths())

