# Native SIAMCAT reference implementation for the fixed source-only benchmark.
# This script is intentionally separate from the Python CLR benchmark because
# SIAMCAT owns its own source-fitted normalization workflow.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 6L) {
  stop(paste(
    "Usage: Rscript run_siamcat_reference.R <source_matrix> <source_metadata>",
    "<target_matrix> <target_metadata_without_labels> <output_dir> <seed>"
  ))
}
source_matrix_path <- arguments[[1]]
source_metadata_path <- arguments[[2]]
target_matrix_path <- arguments[[3]]
target_metadata_path <- arguments[[4]]
output_dir <- arguments[[5]]
seed <- as.integer(arguments[[6]])

if (!requireNamespace("SIAMCAT", quietly = TRUE)) {
  stop("SIAMCAT is not installed. Install it in R_LIBS_USER before running this reference.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

read_matrix <- function(path) {
  data <- read.csv(path, row.names = 1, check.names = FALSE)
  as.matrix(data)
}
read_metadata <- function(path) {
  metadata <- read.csv(path, stringsAsFactors = FALSE)
  rownames(metadata) <- metadata$sample_id
  metadata
}

source_data <- read_matrix(source_matrix_path)
source_meta <- read_metadata(source_metadata_path)
target_data <- read_matrix(target_matrix_path)
target_meta <- read_metadata(target_metadata_path)
if (!"outcome_ibd" %in% colnames(source_meta)) {
  stop("Source metadata must contain outcome_ibd.")
}
if ("outcome_ibd" %in% colnames(target_meta)) {
  stop("Target metadata must not contain outcome_ibd; labels are locked outside R training.")
}
if (!identical(rownames(source_data), rownames(target_data))) {
  stop("Source and target must have identical, ordered shared species.")
}
if (!identical(colnames(source_data), source_meta$sample_id)) {
  stop("Source feature matrix and metadata sample order differ.")
}
if (!identical(colnames(target_data), target_meta$sample_id)) {
  stop("Target feature matrix and metadata sample order differ.")
}

source_labels <- ifelse(source_meta$outcome_ibd == 1L, "IBD", "Control")
names(source_labels) <- source_meta$sample_id
labels <- SIAMCAT::create.label(
  label = source_labels,
  case = "IBD"
)
model <- SIAMCAT::siamcat(feat = source_data, label = labels)
# All filtering and normalization are fitted on the source cohort only.
model <- SIAMCAT::filter.features(model, filter.method = "abundance", cutoff = 1e-06)
model <- SIAMCAT::normalize.features(
  model,
  norm.method = "log.std",
  norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1)
)
# The split is defined only in the source cohort. Five source-trained LASSO
# models are later averaged on the locked external target.
model <- SIAMCAT::create.data.split(model, num.folds = 5, num.resample = 1, stratify = TRUE)
model <- SIAMCAT::train.model(model, method = "lasso")

target_holdout <- SIAMCAT::siamcat(feat = target_data)
prediction <- SIAMCAT::make.predictions(model, siamcat.holdout = target_holdout)
score_matrix <- as.matrix(SIAMCAT::pred_matrix(prediction))
if (nrow(score_matrix) == nrow(target_meta)) {
  scores <- rowMeans(score_matrix)
} else if (ncol(score_matrix) == nrow(target_meta)) {
  scores <- colMeans(score_matrix)
} else {
  stop("SIAMCAT prediction matrix does not align to target samples.")
}
if (length(scores) != nrow(target_meta)) {
  stop("SIAMCAT prediction count does not match target metadata.")
}
result <- data.frame(
  sample_id = target_meta$sample_id,
  probability = scores
)
write.csv(result, file.path(output_dir, "siamcat_external_predictions.csv"), row.names = FALSE)
saveRDS(model, file.path(output_dir, "siamcat_source_model.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "siamcat_session_info.txt"))

