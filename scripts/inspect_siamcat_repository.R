local_library <- Sys.getenv(
  "R_LIBS_USER",
  unset = Sys.getenv("IBD_R_LIB", unset = "")
)
.libPaths(c(local_library, .libPaths()))
options(repos = BiocManager::repositories())
packages <- available.packages()
print(packages["SIAMCAT", c("Package", "Version", "Repository")])

