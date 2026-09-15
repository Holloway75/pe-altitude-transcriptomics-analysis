#!/usr/bin/env Rscript

# MR-JTI runner retaining the vendored basic residual-bootstrap inference.
suppressPackageStartupMessages({
  library(optparse)
  library(glmnet)
  library(HDCI)
})

options_list <- list(
  make_option("--df-path", dest="df_path", type="character"),
  make_option("--vendor-script", dest="vendor_script", type="character"),
  make_option("--result-path", dest="result_path", type="character"),
  make_option("--bootstrap-path", dest="bootstrap_path", type="character"),
  make_option("--n-folds", dest="n_folds", type="integer", default=5L),
  make_option("--n-bootstrap", dest="n_bootstrap", type="integer", default=500L),
  make_option("--n-genes", dest="n_genes", type="integer", default=1L),
  make_option("--min-snps", dest="min_snps", type="integer", default=20L),
  make_option("--seed", type="integer", default=20260714L)
)
opt <- parse_args(OptionParser(option_list=options_list))
if (is.na(opt$n_genes) || opt$n_genes < 1L) stop("--n-genes must be a positive planned family size")

vendor_lines <- readLines(opt$vendor_script, warn=FALSE)
start <- grep("^TRB_LASSO<-", vendor_lines)[1]
end <- grep("^#load df", vendor_lines)[1]
if (is.na(start) || is.na(end) || start >= end) stop("Cannot locate vendored MR-JTI functions")
eval(parse(text=paste(vendor_lines[start:(end - 1L)], collapse="\n")), envir=.GlobalEnv)

df <- read.table(opt$df_path, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
required <- c("rsid", "ldscore", "eqtl_beta", "eqtl_se", "eqtl_p", "gwas_beta", "gwas_se", "gwas_p")
if (!all(required %in% names(df))) stop("MR-JTI input is missing required columns")
if (nrow(df) < opt$min_snps) stop("Need >= min-snps variants to run MR-JTI")
if (any(!is.finite(as.matrix(df[, setdiff(required, "rsid")]))) || anyDuplicated(df$rsid)) stop("Invalid MR-JTI input")

set.seed(opt$seed)
df <- df[order(df$gwas_p, df$rsid), ]
n_gwas <- sum(df$gwas_p < 0.05)
penalty.factor <- c(0, 0, rep(1, n_gwas))
model_df <- df[, c("rsid", "gwas_beta", "eqtl_beta", "ldscore")]
if (n_gwas > 0) model_df[, (ncol(model_df)+1):(ncol(model_df)+n_gwas)] <- diag(nrow(model_df))[, 1:n_gwas, drop=FALSE]
y <- as.numeric(scale(model_df[, "gwas_beta"]))
x <- apply(as.matrix(model_df[, 3:ncol(model_df), drop=FALSE]), 2, scale)
if (any(!is.finite(y)) || any(!is.finite(x))) stop("Zero-variance or nonfinite standardized MR-JTI column")

# The vendored function does not expose bootstrap draws.  This protocol copy
# reproduces it and additionally retains the draws and stability diagnostics.
TRB_with_draws <- function(x, y, B, alpha, nfolds, weights, penalty.factor) {
  x <- as.matrix(x); y <- as.numeric(y); n <- nrow(x); p <- ncol(x)
  globalfit <- glmnet(x, y, standardize=TRUE, intercept=TRUE, weights=weights, penalty.factor=penalty.factor)
  cvfit <- weighted.escv.glmnet(x, y, lambda=globalfit$lambda, nfolds=nfolds,
    tau=0, cv.OLS=FALSE, standardize=TRUE, intercept=TRUE, weights=weights,
    penalty.factor=penalty.factor)
  lambda.opt <- cvfit$lambda.cv
  fit_value <- predict(globalfit, newx=x, s=lambda.opt)
  original_beta <- as.numeric(predict(globalfit, type="coefficients", s=lambda.opt))[-1]
  original_beta <- ifelse(abs(original_beta) < 1/ncol(x), 0, original_beta)
  residual_center <- (y - fit_value) - mean(y - fit_value)
  draws <- matrix(0, nrow=B, ncol=p)
  for (i in seq_len(B)) {
    ystar <- fit_value + residual_center[sample.int(n, n, replace=TRUE)]
    boot <- weighted.Lasso(x=x, y=ystar, lambda=lambda.opt, standardize=TRUE,
      intercept=TRUE, weights=weights, penalty.factor=penalty.factor)
    draws[i, ] <- ifelse(abs(boot$beta) < 1/ncol(x), 0, boot$beta)
  }
  percentile <- apply(draws, 2, quantile, probs=c(1-alpha/2, alpha/2))
  interval <- rbind(2 * original_beta - percentile[1, ], 2 * original_beta - percentile[2, ])
  cv_index <- which.min(abs(cvfit$lambda - lambda.opt))
  list(beta=colMeans(draws), original_beta=original_beta, interval=interval, draws=draws,
    lambda=lambda.opt, cv_mean_error=cvfit$cv[cv_index], cv_error_se=cvfit$cv.error[cv_index])
}

family_alpha <- 0.05 / opt$n_genes
ans <- TRB_with_draws(x, y, opt$n_bootstrap, family_alpha, opt$n_folds, rep(1, nrow(df)), penalty.factor)
draw <- ans$draws[, 1]
ci_significance <- ifelse(ans$interval[1,1] * ans$interval[2,1] > 0, "sig", "nonsig")
result <- data.frame(standardized_expression_coefficient=ans$beta[1],
  original_lasso_coefficient=ans$original_beta[1], ci_low=ans$interval[1,1],
  ci_high=ans$interval[2,1], ci_method="vendored_basic_residual_bootstrap_bonferroni",
  ci_significance=ci_significance, n_genes=opt$n_genes, family_alpha=family_alpha,
  pvalue=NA, qvalue=NA, inference_status="native_bonferroni_ci_selection_dependent",
  n_snps=nrow(df), n_folds=opt$n_folds, n_bootstrap=opt$n_bootstrap,
  seed=opt$seed, lambda=ans$lambda, cv_mean_error=ans$cv_mean_error,
  cv_error_se=ans$cv_error_se, bootstrap_zero_fraction=mean(draw == 0),
  bootstrap_positive_fraction=mean(draw > 0), bootstrap_negative_fraction=mean(draw < 0))
write.table(result, opt$result_path, sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(iteration=seq_along(draw), expression_beta=draw), opt$bootstrap_path,
  sep="\t", quote=FALSE, row.names=FALSE)
