#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(limma)
})

parse_args <- function(x) {
  if (length(x) %% 2L != 0L) stop("Arguments must be --name value pairs")
  setNames(as.list(x[seq(2L, length(x), 2L)]),
           sub("^--", "", x[seq(1L, length(x), 2L)]))
}

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)
args <- parse_args(commandArgs(trailingOnly = TRUE))
project <- normalizePath(args[["project-root"]], mustWork = TRUE)
output <- args[["output-dir"]]
if (!grepl("^/", output)) output <- file.path(project, output)
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty output: ", output)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)

analysis_id <- "m2_gse103927_composition_sensitivity_v3_20260910"
base <- file.path(project, "results/02_altitude_response/analyses",
                  "m2_gse103927_alt16_paired_v2_20260910/GSE103927")
release <- file.path(project, "data/processed/02_altitude_response/releases",
                     "m2_module3_input_freeze_v2_20260910")
v5_locked <- file.path(project, "data/processed/03_bidirectional_convergence/releases",
                       "m3_exploratory_shared_programs_v6_20260910")

sensitivity_path <- file.path(
  release, "composition_diagnostics/GSE103927_zero_proxy_change_gene_sensitivity.tsv"
)
proxy_summary_path <- file.path(base, "cell_composition_proxy_summary.tsv")
proxy_change_path <- file.path(
  release, "composition_diagnostics/GSE103927_cell_proxy_changes.tsv"
)
universe_path <- file.path(v5_locked, "gene_universe_and_mapping.tsv")
signature_registry_path <- file.path(v5_locked, "compact_signature_registry.tsv")
required <- c(sensitivity_path, proxy_summary_path, proxy_change_path,
              universe_path, signature_registry_path)
if (any(!file.exists(required))) {
  stop("Missing required input(s): ", paste(required[!file.exists(required)],
                                             collapse = ", "))
}

x <- fread(sensitivity_path)
proxies <- fread(proxy_summary_path)
universe <- fread(universe_path)
signature_registry <- fread(signature_registry_path)
if (anyDuplicated(x$hgnc_id) || anyDuplicated(universe$hgnc_id)) {
  stop("HGNC IDs must be unique")
}

metric_row <- function(label, z) {
  primary_sig <- z$primary_fdr < 0.05
  conditional_sig <- z$proxy_conditional_fdr < 0.05
  data.table(
    analysis_id = analysis_id,
    gene_domain = label,
    n_genes = nrow(z),
    pearson_log2fc = cor(z$primary_logFC, z$proxy_conditional_logFC,
                         method = "pearson"),
    spearman_log2fc = cor(z$primary_logFC, z$proxy_conditional_logFC,
                          method = "spearman"),
    pearson_t = cor(z$primary_t, z$proxy_conditional_t, method = "pearson"),
    spearman_t = cor(z$primary_t, z$proxy_conditional_t, method = "spearman"),
    sign_concordance = mean(sign(z$primary_logFC) ==
                              sign(z$proxy_conditional_logFC)),
    primary_fdr_lt_0_05 = sum(primary_sig),
    conditional_fdr_lt_0_05 = sum(conditional_sig),
    both_fdr_lt_0_05 = sum(primary_sig & conditional_sig),
    both_fdr_same_direction = sum(primary_sig & conditional_sig &
                                    sign(z$primary_logFC) ==
                                    sign(z$proxy_conditional_logFC)),
    primary_strict_fdr_fc05 = sum(primary_sig & abs(z$primary_logFC) >= 0.5),
    conditional_strict_fdr_fc05 = sum(conditional_sig &
                                        abs(z$proxy_conditional_logFC) >= 0.5)
  )
}

joined <- x[match(universe$hgnc_id, hgnc_id)]
if (anyNA(joined$hgnc_id)) stop("Composition table does not cover v5 universe")
metrics <- rbind(metric_row("all_module2_genes", x),
                 metric_row("v5_8101_gene_universe", joined))
fwrite(metrics, file.path(output, "composition_effect_stability.tsv"),
       sep = "\t", quote = FALSE)

proxy_out <- copy(proxies)
proxy_out[, analysis_id := analysis_id]
setcolorder(proxy_out, c("analysis_id", setdiff(names(proxy_out), "analysis_id")))
fwrite(proxy_out, file.path(output, "composition_proxy_changes_summary.tsv"),
       sep = "\t", quote = FALSE)

strict_primary <- joined[primary_fdr < 0.05 & abs(primary_logFC) >= 0.5]
strict_primary[, `:=`(
  analysis_id = analysis_id,
  conditional_fdr_pass = proxy_conditional_fdr < 0.05,
  conditional_strict_pass = proxy_conditional_fdr < 0.05 &
    abs(proxy_conditional_logFC) >= 0.5,
  direction_preserved = sign(primary_logFC) == sign(proxy_conditional_logFC)
)]
setcolorder(strict_primary, c("analysis_id", setdiff(names(strict_primary),
                                                      "analysis_id")))
fwrite(strict_primary, file.path(output, "primary_strict_gene_sensitivity.tsv"),
       sep = "\t", quote = FALSE)

sets <- lapply(signature_registry$universe_genes,
               function(z) strsplit(z, "|", fixed = TRUE)[[1L]])
names(sets) <- signature_registry$signature_id
indices <- lapply(sets, function(genes) which(joined$gene_symbol %in% genes))
run_camera_pr <- function(statistic, prefix) {
  result <- as.data.table(cameraPR(statistic, indices,
                                    inter.gene.cor = 0.01,
                                    sort = FALSE),
                          keep.rownames = "signature_id")
  setnames(result, c("NGenes", "Direction", "PValue", "FDR"),
           paste0(prefix, c("_n_genes", "_direction", "_p", "_bh_fdr")))
  result
}
primary_signature <- run_camera_pr(joined$primary_t, "primary")
conditional_signature <- run_camera_pr(joined$proxy_conditional_t,
                                       "proxy_conditional")
signature_comparison <- merge(primary_signature, conditional_signature,
                              by = "signature_id", sort = FALSE)
signature_comparison[, `:=`(
  analysis_id = analysis_id,
  direction_preserved = primary_direction == proxy_conditional_direction,
  interpretation = paste0(
    "proxy-conditional intercept estimates expression change at zero change ",
    "in all four marker scores and is an extrapolative sensitivity estimand"
  )
)]
setcolorder(signature_comparison, c("analysis_id", "signature_id",
                                    setdiff(names(signature_comparison),
                                            c("analysis_id", "signature_id"))))
fwrite(signature_comparison,
       file.path(output, "compact_signature_composition_sensitivity.tsv"),
       sep = "\t", quote = FALSE)

wars <- x[gene_symbol == "WARS1"]
if (nrow(wars) != 1L) stop("Expected one WARS1 row")
fwrite(wars, file.path(output, "wars1_composition_sensitivity.tsv"),
       sep = "\t", quote = FALSE)

sources <- data.table(
  source = c("gene_sensitivity", "proxy_summary", "proxy_changes",
             "v5_universe", "compact_signature_registry"),
  path = required,
  bytes = file.info(required)$size,
  sha256 = vapply(required, sha256, character(1L))
)
fwrite(sources, file.path(output, "source_manifest.tsv"), sep = "\t",
       quote = FALSE)

universe_metrics <- metrics[gene_domain == "v5_8101_gene_universe"]
summary <- c(
  "# GSE103927 cell-composition sensitivity summary",
  "",
  paste0("Analysis ID: `", analysis_id, "`"),
  "",
  "The primary estimand remains the total paired PBMC response. The sensitivity",
  "model conditions on changes in B-cell, NK-cell, T-cell and monocyte marker",
  "scores; its intercept is an extrapolative zero-marker-change estimand and is",
  "not a replacement primary analysis.",
  "",
  "## Composition diagnostics",
  "",
  paste0("- T-cell marker score: mean change ",
         sprintf("%.3f", proxies[cell_proxy == "T_cell", mean_change]),
         ", FDR ", sprintf("%.3g", proxies[cell_proxy == "T_cell", wilcoxon_fdr]), "."),
  paste0("- NK-cell marker score: mean change ",
         sprintf("%.3f", proxies[cell_proxy == "NK_cell", mean_change]),
         ", FDR ", sprintf("%.3g", proxies[cell_proxy == "NK_cell", wilcoxon_fdr]), "."),
  paste0("- Monocyte marker score: mean change ",
         sprintf("%.3f", proxies[cell_proxy == "monocyte", mean_change]),
         ", FDR ", sprintf("%.3g", proxies[cell_proxy == "monocyte", wilcoxon_fdr]), "."),
  paste0("- B-cell marker score: mean change ",
         sprintf("%.3f", proxies[cell_proxy == "B_cell", mean_change]),
         ", FDR ", sprintf("%.3g", proxies[cell_proxy == "B_cell", wilcoxon_fdr]), "."),
  "",
  "## Effect stability in the frozen 8,101-gene universe",
  "",
  paste0("- Primary versus conditional log2FC: Pearson r=",
         sprintf("%.3f", universe_metrics$pearson_log2fc),
         ", Spearman rho=", sprintf("%.3f", universe_metrics$spearman_log2fc), "."),
  paste0("- Direction concordance: ",
         sprintf("%.1f%%", 100 * universe_metrics$sign_concordance), "."),
  paste0("- FDR<0.05 genes: primary ", universe_metrics$primary_fdr_lt_0_05,
         "; proxy-conditional ", universe_metrics$conditional_fdr_lt_0_05, "."),
  "",
  "Large marker-score shifts and loss of most gene-level significance after",
  "conditioning show that cell-mixture change is a major component of the total",
  "PBMC response. This does not invalidate the frozen total-PBMC estimand, but it",
  "precludes interpreting all primary DE genes as within-cell regulation."
)
writeLines(summary, file.path(output, "COMPOSITION_SENSITIVITY_SUMMARY.md"))

validation <- list(
  analysis_id = analysis_id,
  status = "PASS",
  checks = list(
    list(id = "unique_gene_rows", pass = !anyDuplicated(x$hgnc_id)),
    list(id = "v5_universe_closed", pass = nrow(joined) == 8101L &&
           !anyNA(joined$hgnc_id)),
    list(id = "strict_gene_table_nonempty", pass = nrow(strict_primary) > 0L),
    list(id = "six_signature_rows", pass = nrow(signature_comparison) == 6L),
    list(id = "four_proxy_rows", pass = nrow(proxies) == 4L),
    list(id = "wars1_unique", pass = nrow(wars) == 1L)
  )
)
if (!all(vapply(validation$checks, `[[`, logical(1L), "pass"))) {
  validation$status <- "FAIL"
}
write_json(validation, file.path(output, "validation_report.json"),
           pretty = TRUE, auto_unbox = TRUE)
if (validation$status != "PASS") stop("Validation failed")
cat(toJSON(list(status = "PASS", output = output), auto_unbox = TRUE), "\n")
