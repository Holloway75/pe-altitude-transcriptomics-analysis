#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

parse_args <- function(x) {
  if (length(x) %% 2L != 0L) stop("Arguments must be --name value pairs")
  setNames(as.list(x[seq(2L, length(x), 2L)]),
           sub("^--", "", x[seq(1L, length(x), 2L)]))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project <- normalizePath(args[["project-root"]], mustWork = TRUE)
output_dir <- args[["output-dir"]]
if (!grepl("^/", output_dir)) output_dir <- file.path(project, output_dir)
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty frozen release: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "donor_contrasts"), showWarnings = FALSE)
dir.create(file.path(output_dir, "gene_statistics"), showWarnings = FALSE)
dir.create(file.path(output_dir, "metadata"), showWarnings = FALSE)
dir.create(file.path(output_dir, "composition_diagnostics"),
           showWarnings = FALSE)

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
                                      na = "NA")

axes <- list(
  list(
    axis_id = "ALT_GSE103927_ALT16_TOTAL",
    dataset_id = "GSE103927",
    analysis_id = "m2_gse103927_alt16_paired_v2_20260910",
    role = "PRIMARY_FORMAL",
    tissue = "PBMC",
    population = "21 healthy lowlanders",
    altitude_m = 5260L,
    exposure_phase = "day_16_acclimatized",
    contrast = "ALT16_5260m_minus_sea_level_pre_ascent_baseline",
    coexposures = "expedition exposure; cell-mixture change",
    object = file.path(
      project, "results/02_altitude_response/analyses",
      "m2_gse103927_alt16_paired_v2_20260910/GSE103927/model_objects.rds"
    ),
    de = file.path(
      project, "results/02_altitude_response/analyses",
      "m2_gse103927_alt16_paired_v2_20260910/GSE103927/differential_expression.tsv"
    ),
    metadata = file.path(
      project, "results/02_altitude_response/analyses",
      "m2_gse103927_alt16_paired_v2_20260910/GSE103927/sample_metadata.tsv"
    ),
    proxy_changes = file.path(
      project, "results/02_altitude_response/analyses",
      "m2_gse103927_alt16_paired_v2_20260910/GSE103927/cell_composition_proxy_changes.tsv"
    ),
    proxy_gene_sensitivity = file.path(
      project, "results/02_altitude_response/analyses",
      "m2_gse103927_alt16_paired_v2_20260910/GSE103927/composition_proxy_sensitivity.tsv"
    ),
    formal = TRUE,
    independent_replication = FALSE
  ),
  list(
    axis_id = "ALT_GSE333506_HAN_WEEK4_TOTAL",
    dataset_id = "GSE333506",
    analysis_id = "m2_gse333506_han_week4_paired_v2_20260910",
    role = "CONTEXT_SENSITIVITY",
    tissue = "whole blood",
    population = "7 Han Chinese male endurance athletes",
    altitude_m = 2260L,
    exposure_phase = "week_4_altitude_training",
    contrast = "week4_2260m_minus_pre_ascent_baseline",
    coexposures = "four weeks endurance training inseparable from altitude",
    object = file.path(
      project, "results/02_altitude_response/analyses",
      "m2_gse333506_han_week4_paired_v2_20260910/GSE333506/model_objects.rds"
    ),
    de = file.path(
      project, "results/02_altitude_response/analyses",
      "m2_gse333506_han_week4_paired_v2_20260910/GSE333506/differential_expression.tsv"
    ),
    metadata = file.path(
      project, "results/02_altitude_response/analyses",
      "m2_gse333506_han_week4_paired_v2_20260910/GSE333506/sample_metadata.tsv"
    ),
    proxy_changes = file.path(
      project, "results/02_altitude_response/analyses",
      "m2_gse333506_han_week4_paired_v2_20260910/GSE333506/cell_composition_proxy_changes.tsv"
    ),
    proxy_gene_sensitivity = NA_character_,
    formal = FALSE,
    independent_replication = FALSE
  )
)

required <- unique(unlist(lapply(axes, function(x) {
  unname(unlist(x[c("object", "de", "metadata", "proxy_changes")]))
})))
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing release source(s): ", paste(missing,
                                                               collapse = ", "))

axis_rows <- list()
source_rows <- list()
validation_rows <- list()

for (axis in axes) {
  model <- readRDS(axis$object)
  difference <- model$paired_difference
  annotation <- as.data.table(model$annotation)
  de <- fread(axis$de)

  if (!is.matrix(difference) || nrow(difference) == 0L ||
      ncol(difference) < 3L) {
    stop(axis$axis_id, ": invalid paired_difference matrix")
  }
  if (anyNA(difference) || any(!is.finite(difference))) {
    stop(axis$axis_id, ": paired_difference contains non-finite values")
  }
  if (anyDuplicated(rownames(difference)) ||
      anyDuplicated(colnames(difference))) {
    stop(axis$axis_id, ": duplicated gene or donor identifiers")
  }
  annotation <- annotation[match(rownames(difference), gene_symbol)]
  if (anyNA(annotation$gene_symbol) || anyNA(annotation$hgnc_id) ||
      anyDuplicated(annotation$hgnc_id)) {
    stop(axis$axis_id, ": HGNC mapping is incomplete or non-unique")
  }
  de <- de[match(rownames(difference), gene_symbol)]
  if (anyNA(de$gene_symbol)) stop(axis$axis_id, ": DE-to-matrix mapping failed")

  observed_mean <- rowMeans(difference)
  max_delta <- max(abs(observed_mean - de$logFC))
  if (!is.finite(max_delta) || max_delta > 1e-8) {
    stop(axis$axis_id, ": donor means do not reproduce module-2 logFC; max ",
         max_delta)
  }

  n <- ncol(difference)
  observed_sd <- apply(difference, 1L, sd)
  observed_se <- observed_sd / sqrt(n)
  observed_t <- observed_mean / observed_se
  observed_t[observed_sd == 0 & observed_mean == 0] <- 0
  observed_t[observed_sd == 0 & observed_mean != 0] <-
    sign(observed_mean[observed_sd == 0 & observed_mean != 0]) * Inf
  probability <- pt(observed_t, df = n - 1L)
  probability <- pmin(1 - .Machine$double.eps,
                      pmax(.Machine$double.xmin, probability))
  signed_z <- qnorm(probability)
  p_two_sided <- 2 * pt(-abs(observed_t), df = n - 1L)

  contrast_table <- as.data.table(difference)
  contrast_table[, hgnc_symbol := annotation$gene_symbol]
  contrast_table[, hgnc_id := annotation$hgnc_id]
  setcolorder(contrast_table,
              c("hgnc_id", "hgnc_symbol",
                setdiff(names(contrast_table), c("hgnc_id", "hgnc_symbol"))))
  setorder(contrast_table, hgnc_id)

  statistics <- data.table(
    axis_id = axis$axis_id,
    hgnc_id = annotation$hgnc_id,
    hgnc_symbol = annotation$gene_symbol,
    n_donors = n,
    mean_difference = observed_mean,
    sd_difference = observed_sd,
    se_difference = observed_se,
    t_one_sample = observed_t,
    signed_z_one_sample = signed_z,
    p_two_sided_one_sample = p_two_sided,
    module2_moderated_logFC = de$logFC,
    module2_moderated_t = de$t,
    module2_p = de$P.Value,
    module2_fdr = de$adj.P.Val
  )
  setorder(statistics, hgnc_id)

  contrast_path <- file.path(output_dir, "donor_contrasts",
                             paste0(axis$axis_id, ".tsv.gz"))
  statistics_path <- file.path(output_dir, "gene_statistics",
                               paste0(axis$axis_id, ".tsv.gz"))
  metadata_path <- file.path(output_dir, "metadata",
                             paste0(axis$dataset_id, "_sample_metadata.tsv"))
  proxy_path <- file.path(
    output_dir, "composition_diagnostics",
    paste0(axis$dataset_id, "_cell_proxy_changes.tsv")
  )
  write_tsv(contrast_table, contrast_path)
  write_tsv(statistics, statistics_path)
  file.copy(axis$metadata, metadata_path, overwrite = FALSE)
  file.copy(axis$proxy_changes, proxy_path, overwrite = FALSE)
  if (!is.na(axis$proxy_gene_sensitivity)) {
    file.copy(
      axis$proxy_gene_sensitivity,
      file.path(output_dir, "composition_diagnostics",
                paste0(axis$dataset_id,
                       "_zero_proxy_change_gene_sensitivity.tsv")),
      overwrite = FALSE
    )
  }

  axis_rows[[length(axis_rows) + 1L]] <- data.table(
    axis_id = axis$axis_id,
    dataset_id = axis$dataset_id,
    analysis_id = axis$analysis_id,
    formal_role = axis$role,
    tissue = axis$tissue,
    population = axis$population,
    altitude_m = axis$altitude_m,
    exposure_phase = axis$exposure_phase,
    contrast = axis$contrast,
    positive_direction = "higher expression after/during altitude exposure",
    n_donors = n,
    n_hgnc_genes = nrow(difference),
    paired_design = TRUE,
    coexposures = axis$coexposures,
    composition_status = if (axis$formal)
      "total_PBMC_primary; proxy_change_diagnostic_only" else
      "total_whole_blood_context; proxy_change_diagnostic_only",
    eligible_for_formal_alt_axis = axis$formal,
    eligible_as_independent_replication = axis$independent_replication,
    donor_contrast_path = sub(paste0("^", output_dir, "/"), "",
                              contrast_path),
    gene_statistics_path = sub(paste0("^", output_dir, "/"), "",
                               statistics_path)
  )

  for (source in c("object", "de", "metadata", "proxy_changes")) {
    path <- axis[[source]]
    source_rows[[length(source_rows) + 1L]] <- data.table(
      dataset_id = axis$dataset_id,
      source_type = source,
      path = normalizePath(path),
      bytes = file.info(path)$size,
      sha256 = sha256(path)
    )
  }
  validation_rows[[length(validation_rows) + 1L]] <- data.table(
    axis_id = axis$axis_id,
    check_id = c("finite_donor_matrix", "unique_hgnc_id",
                 "unique_donor_id", "mean_reproduces_module2_logFC",
                 "positive_direction_frozen"),
    status = "PASS",
    detail = c(
      sprintf("%d genes x %d donors", nrow(difference), n),
      sprintf("%d unique HGNC IDs", nrow(difference)),
      sprintf("%d unique donors", n),
      sprintf("maximum absolute difference %.3g", max_delta),
      "all contrasts are altitude/post-exposure minus pre/low-altitude"
    )
  )
}

axis_registry <- rbindlist(axis_rows)
source_manifest <- rbindlist(source_rows)
validation <- rbindlist(validation_rows)
write_tsv(axis_registry, file.path(output_dir, "axis_registry.tsv"))
write_tsv(source_manifest, file.path(output_dir, "source_manifest.tsv"))
write_tsv(validation, file.path(output_dir, "validation_report.tsv"))

excluded <- data.table(
  dataset_id = c("GSE196728", "GSE103940", "GSE46480", "GSE100988",
                 "GSE75665"),
  module3_status = "EXCLUDED_ARCHIVED",
  reason = c(
    paste("RI 5100 m is contrasted with a post-expedition sea-level visit;",
          "fixed visit order, medication and expedition effects are inseparable"),
    paste("acute whole-blood arrival response; current input is submitted FPKM;",
          "not the acclimatized pre-ascent-baseline estimand"),
    paste("South Pole day-3 environmental response with ambiguous blood/PBMC",
          "material and unresolvable travel, cold, circadian and medication effects"),
    paste("independent term placenta samples estimate long-term residence/ancestry",
          "context, not an adult paired blood response"),
    "retired previously; acute/AMS study and no longer in the active result set"
  )
)
write_tsv(excluded, file.path(output_dir, "excluded_datasets.tsv"))

registry <- list(
  schema_version = "1.0",
  release_id = "m2_module3_input_freeze_v2_20260910",
  status = "FROZEN_MODULE2_INPUT_FOR_MODULE3_DEVELOPMENT",
  frozen_at = "2026-09-10",
  revision_note = paste(
    "v2 (2026-09-10): module-2 eBayes trend covariate corrected to mean log2",
    "expression abundance (v1 regressed the prior variance on the row mean of the",
    "difference matrix, equal to logFC under the intercept design); see",
    "docs/report/PRE_SUBMISSION_CODE_AUDIT_20260910.md (finding A-1).",
    "Donor contrast matrices are numerically unchanged; module2_moderated_*",
    "columns are recomputed under the corrected model."
  ),
  research_estimand =
    "within-person acclimatized altitude response relative to pre-ascent baseline",
  formal_alt_axis = "ALT_GSE103927_ALT16_TOTAL",
  context_sensitivity_axis = "ALT_GSE333506_HAN_WEEK4_TOTAL",
  independent_replication_count = 0,
  formal_inference_boundary = paste(
    "GSE333506 is not an independent confirmatory replicate and must not enter",
    "ALT partial-conjunction or formal joint P values."
  ),
  composition_boundary = paste(
    "Total PBMC is the formal GSE103927 estimand. Marker-score changes and the",
    "zero-proxy-change intercept are diagnostics only and are not independent axes."
  ),
  gene_universe_rule = paste(
    "Each axis retains all unique current-HGNC genes available in that study.",
    "The formal PE-JTI intersection must be built downstream from GSE103927 only;",
    "context datasets must not shrink the formal universe."
  ),
  axis_registry = "axis_registry.tsv",
  excluded_datasets = "excluded_datasets.tsv",
  source_manifest = "source_manifest.tsv",
  validation_report = "validation_report.tsv"
)
write_json(registry, file.path(output_dir, "freeze_registry.json"),
           pretty = TRUE, auto_unbox = TRUE)

readme <- c(
  "# Frozen Module 2 input for Module 3",
  "",
  "Release: `m2_module3_input_freeze_v2_20260910`",
  "",
  "GSE103927 ALT16 minus pre-ascent sea-level baseline is the only formal ALT",
  "axis. GSE333506 week 4 minus pre-ascent baseline is retained only as an",
  "external adapted-layer sensitivity because altitude and endurance training",
  "are inseparable and n=7.",
  "",
  "The donor-contrast matrices preserve the experimental unit required by the",
  "Module 3 continuous donor-space rotation. Positive values always mean higher",
  "expression after/during altitude exposure.",
  "",
  "GSE333506 must not be counted as an independent confirmatory replication.",
  "Composition-proxy files are diagnostic and must not be treated as independent",
  "ALT studies. Context datasets listed in `excluded_datasets.tsv` must not",
  "shrink the formal GSE103927/PE-JTI gene universe."
)
writeLines(readme, file.path(output_dir, "FREEZE_README.md"))

generated <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
generated <- generated[file.info(generated)$isdir %in% FALSE]
manifest <- data.table(
  relative_path = sub(paste0("^", output_dir, "/"), "", generated),
  bytes = file.info(generated)$size,
  sha256 = vapply(generated, sha256, character(1L))
)
setorder(manifest, relative_path)
write_tsv(manifest, file.path(output_dir, "MANIFEST.sha256.tsv"))

message("FREEZE PASS: ", nrow(axis_registry), " axes; formal=",
        registry$formal_alt_axis, "; release=", output_dir)
