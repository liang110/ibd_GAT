options(warn = 1)

if (!requireNamespace("SIAMCAT", quietly = TRUE)) {
  stop("SIAMCAT is not installed")
}

library(SIAMCAT)
description <- utils::packageDescription("SIAMCAT")
installed_db <- utils::installed.packages()
installed <- rownames(installed_db)
dependencies <- unique(unlist(tools::package_dependencies(
  "SIAMCAT",
  db = installed_db,
  which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE
)))
missing <- setdiff(dependencies, installed)

cat("R version: ", R.version.string, "\n", sep = "")
cat("SIAMCAT version: ", description$Version, "\n", sep = "")
cat("SIAMCAT library: ", find.package("SIAMCAT"), "\n", sep = "")
cat("Runtime dependency count: ", length(dependencies), "\n", sep = "")
cat("Missing runtime dependencies: ", length(missing), "\n", sep = "")
if (length(missing) > 0L) {
  writeLines(sort(missing))
  quit(status = 2L)
}
cat("SIAMCAT environment audit: PASS\n")
print(utils::sessionInfo())

