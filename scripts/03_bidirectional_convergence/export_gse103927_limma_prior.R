#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: export_gse103927_limma_prior.R MODEL_RDS OUTPUT_TSV")
model_path <- normalizePath(args[[1L]], mustWork = TRUE)
output_path <- args[[2L]]
if (file.exists(output_path)) stop("Refusing to overwrite: ", output_path)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

model <- readRDS(model_path)
fit <- model$primary_fit
annotation <- as.data.table(model$annotation)
if (is.null(fit$df.prior) || is.null(fit$s2.prior) ||
    is.null(fit$df.residual) || is.null(fit$stdev.unscaled)) {
  stop("Primary fit is missing empirical-Bayes fields")
}
if (!identical(rownames(fit$coefficients), annotation$gene_symbol)) {
  stop("Fit and annotation row order differ")
}

out <- data.table(
  hgnc_id = annotation$hgnc_id,
  hgnc_symbol = annotation$gene_symbol,
  df_residual = as.numeric(fit$df.residual),
  df_prior = as.numeric(fit$df.prior),
  s2_prior = as.numeric(fit$s2.prior),
  stdev_unscaled = as.numeric(fit$stdev.unscaled[, 1L]),
  observed_sigma2 = as.numeric(fit$sigma)^2,
  observed_s2_post = as.numeric(fit$s2.post),
  observed_moderated_t = as.numeric(fit$t[, 1L])
)
if (anyDuplicated(out$hgnc_id) || any(!is.finite(as.matrix(out[, -c(1, 2)])))) {
  stop("Invalid exported prior table")
}
fwrite(out, output_path, sep = "\t", quote = FALSE)
cat("LIMMA PRIOR EXPORT PASS: ", nrow(out), " genes\n", sep = "")

