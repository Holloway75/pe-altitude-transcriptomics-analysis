#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(limma)
})

parse_args <- function(x) {
  if (length(x) %% 2L || any(!startsWith(x[seq(1L, length(x), 2L)], "--"))) {
    stop("Arguments must be --name value pairs", call. = FALSE)
  }
  setNames(as.list(x[seq(2L, length(x), 2L)]),
           sub("^--", "", x[seq(1L, length(x), 2L)]))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("geo-dir", "hgnc", "frozen-dir", "output-dir")
missing <- setdiff(required, names(args))
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

dataset <- "GSE333506"
analysis_id <- "m2_gse333506_han_week4_paired_v2_20260910"
geo_dir <- normalizePath(args[["geo-dir"]], mustWork = TRUE)
hgnc_path <- normalizePath(args[["hgnc"]], mustWork = TRUE)
frozen_dir <- normalizePath(args[["frozen-dir"]], mustWork = FALSE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
dir.create(frozen_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(geo_dir, "GSE333506_series_matrix.txt.gz")
raw_tar <- file.path(geo_dir, "GSE333506_RAW.tar")
scan_zip <- file.path(geo_dir, "GSE333506_Raw_Image_Scan.zip")
stopifnot(file.exists(matrix_path), file.exists(raw_tar), file.exists(scan_zip))

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

sha256 <- function(path) {
  z <- system2("sha256sum", path, stdout = TRUE)
  sub("  .*", "", z[[1L]])
}

quoted_values <- function(line) {
  z <- regmatches(line, gregexpr('"[^"]*"', line, perl = TRUE))[[1L]]
  sub('"$', "", sub('^"', "", z))
}

read_series_metadata <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  lines <- character()
  repeat {
    z <- readLines(con, n = 1000L, warn = FALSE)
    if (!length(z)) break
    lines <- c(lines, z)
    if (any(startsWith(z, "!series_matrix_table_begin"))) break
  }
  get_one <- function(prefix) {
    z <- lines[startsWith(lines, prefix)]
    if (length(z) != 1L) stop("Expected one line for ", prefix)
    quoted_values(z)
  }
  gsm <- get_one("!Sample_geo_accession")
  meta <- data.table(
    sample_id = gsm,
    title = get_one("!Sample_title"),
    description = get_one("!Sample_description")
  )
  characteristics <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
  for (line in characteristics) {
    values <- quoted_values(line)
    if (length(values) != length(gsm)) stop("Characteristic length mismatch")
    key <- sub(":.*$", "", values[[1L]])
    val <- trimws(sub("^[^:]+:", "", values))
    if (key %in% names(meta)) stop("Duplicate characteristic key: ", key)
    meta[, (key) := val]
  }
  required_fields <- c("subject_id", "population", "time point", "gender", "tissue")
  if (length(setdiff(required_fields, names(meta)))) {
    stop("Missing sample characteristics: ",
         paste(setdiff(required_fields, names(meta)), collapse = ", "))
  }
  setnames(meta, "time point", "timepoint")
  meta
}

read_expression <- function(path) {
  x <- fread(path, skip = "!series_matrix_table_begin", header = TRUE,
             showProgress = FALSE)
  if (ncol(x) != 43L) stop("Expected ID_REF plus 42 expression columns")
  feature_id <- as.character(x[[1L]])
  x[[1L]] <- NULL
  expr <- as.matrix(x)
  storage.mode(expr) <- "double"
  rownames(expr) <- feature_id
  if (any(!is.finite(expr))) stop("Non-finite expression values")
  expr
}

read_platform <- function(path) {
  td <- tempfile("gse333506-platform-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  target <- "GPL19615_probe_sequences.txt.gz"
  utils::untar(path, files = target, exdir = td)
  p <- fread(file.path(td, target),
             skip = "ID_REF\tRow\tCol\tControlType", header = TRUE,
             fill = TRUE, showProgress = FALSE)
  required_fields <- c("ControlType", "Probe ID", "Transcript")
  if (length(setdiff(required_fields, names(p)))) {
    stop("GPL19615 fields missing: ",
         paste(setdiff(required_fields, names(p)), collapse = ", "))
  }
  p <- p[, .(
    feature_id = as.character(`Probe ID`),
    control_type = as.integer(ControlType),
    transcript_raw = as.character(Transcript)
  )]
  p <- p[feature_id != "" & !is.na(feature_id)]
  p[, transcript_accession := sub("\\..*$", "", transcript_raw)]
  # Replicate spots share the same Probe ID. Their annotation must agree.
  conflicts <- p[, .(n_control = uniqueN(control_type),
                     n_transcript = uniqueN(transcript_accession)),
                 by = feature_id][n_control > 1L | n_transcript > 1L]
  if (nrow(conflicts)) stop("Conflicting GPL annotation for replicated probes")
  unique(p, by = "feature_id")
}

read_raw_array_qc <- function(path, selected_sample_ids) {
  members <- utils::untar(path, list = TRUE)
  wanted <- members[grepl(
    paste0("^(", paste(selected_sample_ids, collapse = "|"), ")_"),
    basename(members)
  )]
  if (length(wanted) != length(selected_sample_ids)) {
    stop("Could not identify one raw feature file per primary sample")
  }
  td <- tempfile("gse333506-raw-qc-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::untar(path, files = wanted, exdir = td)
  rbindlist(lapply(file.path(td, wanted), function(f) {
    con <- gzfile(f, "rt")
    lines <- readLines(con, n = 20L, warn = FALSE)
    close(con)
    stat_header_i <- which(startsWith(lines, "STATS\t"))[[1L]]
    stat_data_i <- which(seq_along(lines) > stat_header_i &
                           startsWith(lines, "DATA\t"))[[1L]]
    fields <- strsplit(lines[[stat_header_i]], "\t", fixed = TRUE)[[1L]][-1L]
    values <- strsplit(lines[[stat_data_i]], "\t", fixed = TRUE)[[1L]][-1L]
    names(values) <- fields
    getv <- function(name) if (name %in% names(values)) values[[name]] else NA_character_
    in_range <- values[grepl("_IsInRange$", names(values))]
    data.table(
      sample_id = sub("_.*$", "", basename(f)),
      raw_file = basename(f),
      is_good_grid = getv("IsGoodGrid"),
      extraction_status = getv("ExtractionStatus"),
      qc_metric_results = getv("QCMetricResults"),
      noncontrol_well_above_background =
        as.numeric(getv("gNonCtrlNumWellAboveBG")),
      percent_feature_nonuniform =
        as.numeric(getv("AnyColorPrcntFeatNonUnifOL")),
      metrics_out_of_range = sum(in_range != "1", na.rm = TRUE)
    )
  }), use.names = TRUE)
}

build_hgnc_refseq_map <- function(path) {
  h <- fread(path, select = c("hgnc_id", "symbol", "name", "status",
                              "locus_group", "refseq_accession"),
             showProgress = FALSE)
  h <- h[status == "Approved" & !is.na(refseq_accession) &
           refseq_accession != ""]
  h[, transcript_accession := sub("\\..*$", "", refseq_accession)]
  ambiguity <- h[, .(n_symbols = uniqueN(symbol)), by = transcript_accession]
  valid <- ambiguity[n_symbols == 1L, transcript_accession]
  h <- h[transcript_accession %in% valid]
  unique(h[, .(transcript_accession, hgnc_id, gene_symbol = symbol,
               gene_description = name, locus_group)],
         by = "transcript_accession")
}

meta <- read_series_metadata(matrix_path)
expr_all <- read_expression(matrix_path)
if (!identical(colnames(expr_all), meta$sample_id)) {
  stop("Series-matrix sample order does not match metadata")
}

meta[, include_primary := population == "Han Chinese" &
       timepoint %in% c("baseline", "week 4")]
meta[, exclusion_reason := fifelse(
  include_primary, NA_character_,
  fifelse(population != "Han Chinese", "Tibetan population excluded",
          "week 1 excluded from primary contrast")
)]
primary_meta <- meta[include_primary == TRUE]
pair_check <- primary_meta[, .(
  n = .N,
  n_baseline = sum(timepoint == "baseline"),
  n_week4 = sum(timepoint == "week 4")
), by = subject_id]
if (nrow(pair_check) != 7L ||
    any(pair_check$n != 2L | pair_check$n_baseline != 1L |
        pair_check$n_week4 != 1L)) {
  stop("Primary pairing contract failed")
}

baseline <- primary_meta[timepoint == "baseline"][order(subject_id)]
week4 <- primary_meta[timepoint == "week 4"][order(subject_id)]
if (!identical(baseline$subject_id, week4$subject_id)) stop("Pair order mismatch")

platform <- read_platform(raw_tar)
hgnc_map <- build_hgnc_refseq_map(hgnc_path)
ann <- merge(
  data.table(feature_id = rownames(expr_all)),
  platform, by = "feature_id", all.x = TRUE, sort = FALSE
)
ann <- merge(ann, hgnc_map, by = "transcript_accession",
             all.x = TRUE, sort = FALSE)
ann <- ann[match(rownames(expr_all), feature_id)]
ann[, mapping_status := fifelse(
  is.na(control_type), "probe_absent_from_GPL_table",
  fifelse(control_type != 0L, "control_probe",
    fifelse(is.na(transcript_accession) | transcript_accession == "",
            "no_transcript_accession",
      fifelse(is.na(gene_symbol), "refseq_not_in_current_HGNC",
              "mapped_current_HGNC")))
)]

selected_ids <- c(baseline$sample_id, week4$sample_id)
raw_array_qc <- read_raw_array_qc(raw_tar, selected_ids)
expr_selected_probe <- expr_all[, selected_ids, drop = FALSE]
expressed_probe <- rowSums(expr_selected_probe >= 5) >=
  ceiling(ncol(expr_selected_probe) / 2)
analysis_probe <- ann$mapping_status == "mapped_current_HGNC" & expressed_probe
analysis_probe[is.na(analysis_probe)] <- FALSE
ann[, expressed_primary := expressed_probe]
ann[, retained_for_gene_collapse := analysis_probe]
if (sum(analysis_probe) < 5000L) stop("Too few mapped, expressed probes")

expr_probe <- expr_all[analysis_probe, primary_meta$sample_id, drop = FALSE]
ann_probe <- ann[analysis_probe]
expr_gene <- avereps(expr_probe, ID = ann_probe$gene_symbol)
gene_ann <- ann_probe[, .(
  hgnc_id = hgnc_id[[1L]],
  gene_description = gene_description[[1L]],
  locus_group = locus_group[[1L]],
  n_retained_probes = .N,
  retained_probe_ids = paste(feature_id, collapse = "|")
), by = gene_symbol]
gene_ann <- gene_ann[match(rownames(expr_gene), gene_symbol)]

diff <- expr_gene[, week4$sample_id, drop = FALSE] -
  expr_gene[, baseline$sample_id, drop = FALSE]
colnames(diff) <- baseline$subject_id
# 2026-09-10 fix: the trend covariate must be gene-expression abundance
# (mean log2 intensity across all 14 primary-analysis samples), not the row
# mean of the difference matrix (under an intercept-only design the latter
# equals logFC itself; see docs/report/PRE_SUBMISSION_CODE_AUDIT_20260910.md,
# finding A-1).
stopifnot(identical(rownames(diff), rownames(expr_gene)))
abundance <- unname(rowMeans(expr_gene))
design <- matrix(1, nrow = ncol(diff), ncol = 1,
                 dimnames = list(colnames(diff), "week4_minus_baseline"))
fit <- eBayes(lmFit(diff, design), trend = abundance, robust = TRUE)
tt <- as.data.table(topTable(fit, coef = 1, number = Inf, sort.by = "none"),
                    keep.rownames = "gene_symbol")
tt[, AveExpr := NULL]
tt <- merge(gene_ann, tt, by = "gene_symbol", all.y = TRUE, sort = FALSE)
tt[, mean_expression := rowMeans(expr_gene)[match(gene_symbol,
                                                  rownames(expr_gene))]]
tt[, direction := fifelse(logFC > 0, "up", "down")]
tt[, significant_fdr_0.05 := adj.P.Val < 0.05]
tt[, significant_fdr_fc_1.5 :=
     adj.P.Val < 0.05 & abs(logFC) >= log2(1.5)]
setcolorder(tt, c("gene_symbol", "hgnc_id", "gene_description", "locus_group",
                  "n_retained_probes", "retained_probe_ids", "logFC",
                  "mean_expression", "t", "P.Value", "adj.P.Val", "B",
                  "direction", "significant_fdr_0.05",
                  "significant_fdr_fc_1.5"))
setorder(tt, P.Value)

# Leave-one-donor-out sensitivity is descriptive and does not replace the
# seven-pair primary model.
loo <- rbindlist(lapply(seq_len(ncol(diff)), function(i) {
  f <- eBayes(lmFit(diff[, -i, drop = FALSE],
                    matrix(1, nrow = ncol(diff) - 1L, ncol = 1L)),
              trend = abundance, robust = TRUE)
  z <- as.data.table(topTable(f, coef = 1, number = Inf, sort.by = "none"),
                     keep.rownames = "gene_symbol")
  z[, omitted_donor := colnames(diff)[i]]
  z
}))
loo_summary <- loo[, .(
  loo_min_logFC = min(logFC),
  loo_max_logFC = max(logFC),
  loo_sign_concordance = mean(sign(logFC) ==
                                sign(tt$logFC[match(gene_symbol,
                                                    tt$gene_symbol)])),
  loo_max_p = max(P.Value),
  loo_min_p = min(P.Value)
), by = gene_symbol]
loo_summary <- merge(
  tt[, .(gene_symbol, primary_logFC = logFC, primary_p = P.Value,
         primary_fdr = adj.P.Val)],
  loo_summary, by = "gene_symbol", sort = FALSE
)
setorder(loo_summary, primary_p)

marker_sets <- list(
  T_cell = c("CD3D", "CD3E", "TRAC", "IL7R"),
  B_cell = c("CD79A", "MS4A1", "CD37", "CD22"),
  monocyte = c("LYZ", "CTSD", "FCGR3A", "LST1"),
  NK_cell = c("NKG7", "GNLY", "KLRD1", "PRF1"),
  neutrophil = c("FCGR3B", "CSF3R", "S100A8", "S100A9"),
  erythroid = c("HBA1", "HBA2", "HBB", "ALAS2"),
  platelet = c("PPBP", "PF4", "NRGN", "RGS18")
)
score_list <- lapply(names(marker_sets), function(s) {
  k <- intersect(marker_sets[[s]], rownames(expr_gene))
  if (length(k) < 2L) return(NULL)
  data.table(sample_id = colnames(expr_gene), cell_proxy = s,
             score = colMeans(expr_gene[k, , drop = FALSE]),
             marker_count = length(k))
})
scores <- rbindlist(score_list, use.names = TRUE)
scores <- merge(scores, primary_meta[, .(sample_id, subject_id, timepoint)],
                by = "sample_id")
score_wide <- dcast(scores, subject_id + timepoint ~ cell_proxy,
                    value.var = "score")
score_base <- score_wide[timepoint == "baseline"][order(subject_id)]
score_week4 <- score_wide[timepoint == "week 4"][order(subject_id)]
proxy_names <- setdiff(names(score_wide), c("subject_id", "timepoint"))
proxy_change <- as.data.table(
  as.matrix(score_week4[, ..proxy_names] - score_base[, ..proxy_names])
)
proxy_change[, subject_id := score_base$subject_id]
setcolorder(proxy_change, c("subject_id", proxy_names))

sample_out <- meta[, .(
  dataset_id = dataset,
  sample_id,
  expression_column = sample_id,
  donor_id = subject_id,
  population,
  timepoint,
  include_primary,
  exclusion_reason,
  gender,
  tissue,
  assay = "GPL19615 CapitalBio/Agilent lncRNA+mRNA V4 one-color array",
  altitude_m = fifelse(population == "Han Chinese" & timepoint != "baseline",
                       2260L,
                       fifelse(population == "Tibetan" & timepoint == "baseline",
                               3500L, NA_integer_)),
  visit = timepoint,
  pair_id = subject_id,
  source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333506",
  geo_title = title,
  description
)]
pair_out <- data.table(
  donor_id = baseline$subject_id,
  baseline_sample = baseline$sample_id,
  baseline_title = baseline$title,
  week4_sample = week4$sample_id,
  week4_title = week4$title
)

annotation_audit <- ann[, .(
  n_features = .N,
  n_expressed_primary = sum(expressed_primary),
  n_retained_for_gene_collapse = sum(retained_for_gene_collapse)
), by = mapping_status]
setorder(annotation_audit, mapping_status)

source_paths <- c(matrix_path, raw_tar, scan_zip, hgnc_path)
source_manifest <- data.table(
  dataset_id = dataset,
  path = source_paths,
  bytes = file.info(source_paths)$size,
  sha256 = vapply(source_paths, sha256, character(1L))
)

freeze <- list(
  analysis_id = analysis_id,
  dataset_id = dataset,
  status = "FREEZE_RESULT",
  frozen_at = "2026-09-10",
  model_revision = paste(
    "v2 (2026-09-10): eBayes trend covariate corrected from the paired-difference",
    "row mean (equal to logFC under the intercept design) to mean log2 expression",
    "abundance across the 14 primary samples; see",
    "docs/report/PRE_SUBMISSION_CODE_AUDIT_20260910.md (finding A-1).",
    "Donor difference matrices, logFC and ordinary one-sample statistics are",
    "unchanged relative to v1; moderated t/P/FDR shift modestly."
  ),
  population = "Han Chinese endurance athletes",
  n_samples_total = nrow(meta),
  n_samples_primary = nrow(primary_meta),
  n_donors_primary = nrow(pair_check),
  excluded_population = "Tibetan",
  excluded_timepoint = "week 1",
  primary_contrast = "week 4 at approximately 2260 m - pre-ascent baseline",
  direction = "positive logFC means higher expression at week 4",
  value_type = "GEO background-corrected, quantile-normalized log2 intensity",
  feature_mapping =
    "GPL19615 RefSeq transcript accession to approved HGNC refseq_accession",
  independent_filter =
    "mapped non-control probes with intensity >= 5 in at least 7 of 14 primary samples",
  probe_collapse = "mean normalized intensity across retained probes per HGNC symbol",
  primary_model = paste(
    "one-sample empirical-Bayes limma (robust, trend on mean log2 expression",
    "abundance across the 14 primary samples) on seven within-donor",
    "week4-baseline differences"
  ),
  multiple_testing = "Benjamini-Hochberg",
  significance_threshold = "FDR < 0.05",
  composition_rule =
    "marker-score diagnostics only; no covariate adjustment because n=7",
  confirmatory_boundary =
    "training and altitude are inseparable; no AMS scores; 2260 m is lower than GSE196728"
)

write_tsv(sample_out, file.path(frozen_dir, "sample_metadata.tsv"))
write_tsv(pair_out, file.path(frozen_dir, "donor_altitude_mapping.tsv"))
write_tsv(annotation_audit, file.path(frozen_dir, "probe_mapping_audit.tsv"))
write_tsv(gene_ann, file.path(frozen_dir, "gene_annotation.tsv"))
write_tsv(source_manifest, file.path(frozen_dir, "source_manifest.tsv"))
write_json(freeze, file.path(frozen_dir, "FREEZE_RESULT.json"),
           pretty = TRUE, auto_unbox = TRUE)

write_tsv(tt, file.path(output_dir, "differential_expression.tsv"))
write_tsv(tt[1:min(100L, .N)], file.path(output_dir, "top_genes.tsv"))
write_tsv(loo_summary, file.path(output_dir, "leave_one_donor_out.tsv"))
write_tsv(scores, file.path(output_dir, "cell_composition_proxy_scores.tsv"))
write_tsv(proxy_change,
          file.path(output_dir, "cell_composition_proxy_changes.tsv"))
write_tsv(sample_out, file.path(output_dir, "sample_metadata.tsv"))
write_tsv(pair_out, file.path(output_dir, "donor_altitude_mapping.tsv"))
write_tsv(annotation_audit, file.path(output_dir, "probe_mapping_audit.tsv"))
write_tsv(source_manifest, file.path(output_dir, "source_manifest.tsv"))
write_json(c(freeze, list(analysis_status = "complete")),
           file.path(output_dir, "analysis_registry.json"),
           pretty = TRUE, auto_unbox = TRUE)

summary_dt <- data.table(
  dataset_id = dataset,
  contrast = "han_week4_minus_pre_ascent_baseline",
  n_pairs = ncol(diff),
  n_genes_tested = nrow(tt),
  fdr_0.05 = sum(tt$adj.P.Val < 0.05),
  up_fdr_0.05 = sum(tt$adj.P.Val < 0.05 & tt$logFC > 0),
  down_fdr_0.05 = sum(tt$adj.P.Val < 0.05 & tt$logFC < 0),
  fdr_fc_1.5 = sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= log2(1.5)),
  nominal_p_0.05 = sum(tt$P.Value < 0.05)
)
write_tsv(summary_dt, file.path(output_dir, "contrast_summary.tsv"))

saveRDS(
  list(
    analysis_id = analysis_id,
    metadata = primary_meta,
    annotation = gene_ann,
    normalized_expression = expr_gene,
    paired_difference = diff,
    primary_design = design,
    primary_fit = fit,
    cell_proxy_changes = proxy_change
  ),
  file.path(output_dir, "model_objects.rds"),
  compress = "xz"
)

pca <- prcomp(t(expr_gene), center = TRUE, scale. = FALSE)
pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_dt <- data.table(sample_id = rownames(pca$x), PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2])
pca_dt <- merge(pca_dt,
                primary_meta[, .(sample_id, subject_id, timepoint)],
                by = "sample_id")
write_tsv(pca_dt, file.path(output_dir, "pca_scores.tsv"))

pair_cor <- vapply(seq_len(nrow(baseline)), function(i) {
  cor(expr_gene[, baseline$sample_id[i]], expr_gene[, week4$sample_id[i]])
}, numeric(1L))
pair_cor_dt <- rbind(
  data.table(sample_id = baseline$sample_id,
             within_pair_correlation = pair_cor),
  data.table(sample_id = week4$sample_id,
             within_pair_correlation = pair_cor)
)
sample_cor <- cor(expr_gene)
median_profile_cor <- apply(sample_cor, 2L, median)
npc_qc <- min(5L, ncol(pca$x))
pca_scaled <- scale(pca$x[, seq_len(npc_qc), drop = FALSE])
pca_distance <- sqrt(rowSums(pca_scaled^2))
pca_cutoff <- median(pca_distance) + 3 * mad(pca_distance)
sample_qc <- data.table(
  sample_id = colnames(expr_gene),
  expression_median = apply(expr_gene, 2L, median),
  expression_iqr = apply(expr_gene, 2L, IQR),
  median_correlation_to_other_samples = median_profile_cor,
  pca_scaled_distance_first5 = pca_distance,
  pca_outlier_first5 = pca_distance > pca_cutoff
)
sample_qc <- merge(sample_qc, pair_cor_dt, by = "sample_id")
sample_qc <- merge(sample_qc, raw_array_qc, by = "sample_id", all.x = TRUE)
sample_qc <- merge(
  sample_qc,
  primary_meta[, .(sample_id, subject_id, timepoint)],
  by = "sample_id"
)
write_tsv(sample_qc, file.path(output_dir, "sample_qc.tsv"))
p <- ggplot(pca_dt, aes(PC1, PC2, color = timepoint, group = subject_id)) +
  geom_line(alpha = 0.35, color = "grey60") +
  geom_point(size = 2.5) +
  labs(x = sprintf("PC1 (%.1f%%)", pct[1]),
       y = sprintf("PC2 (%.1f%%)", pct[2]),
       title = "GSE333506 Han-athlete primary-analysis sample PCA") +
  theme_bw(base_size = 11)
ggsave(file.path(output_dir, "pca.png"), p, width = 7, height = 5, dpi = 160)

box_dt <- as.data.table(expr_gene, keep.rownames = "gene_symbol")
box_dt <- melt(box_dt, id.vars = "gene_symbol",
               variable.name = "sample_id", value.name = "expression")
box_dt <- merge(box_dt,
                primary_meta[, .(sample_id, subject_id, timepoint)],
                by = "sample_id")
pb <- ggplot(box_dt, aes(sample_id, expression, fill = timepoint)) +
  geom_boxplot(outlier.size = 0.2) +
  labs(x = NULL, y = "GEO normalized log2 intensity",
       title = "GSE333506 primary-analysis expression distribution") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(output_dir, "normalized_expression_boxplot.png"),
       pb, width = 9, height = 5, dpi = 160)

volcano <- copy(tt)
volcano[, neglog10p := -log10(pmax(P.Value, .Machine$double.xmin))]
pv <- ggplot(volcano, aes(logFC, neglog10p,
                          color = significant_fdr_0.05)) +
  geom_point(alpha = 0.55, size = 1) +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#777777")) +
  labs(x = "Paired log2 FC (week 4 - baseline)", y = "-log10(P)",
       title = "GSE333506 Han-athlete paired differential expression") +
  theme_bw(base_size = 11) + guides(color = "none")
ggsave(file.path(output_dir, "volcano.png"), pv, width = 7, height = 5,
       dpi = 160)

top5 <- tt[1:min(5L, .N),
           paste0(gene_symbol, " (logFC=", sprintf("%.3f", logFC),
                  ", FDR=", format(adj.P.Val, digits = 3), ")")]
summary_lines <- c(
  "# GSE333506 Han-athlete week-4 paired differential-expression summary",
  "",
  "## Primary estimand",
  "",
  "Included only the 7 Han Chinese male endurance athletes. The primary comparison is training week 4 at about 2 260 m minus the pre-altitude low-altitude baseline within the same subject; week-1 samples and all Tibetan-athlete samples do not enter the primary model.",
  "",
  "## Results",
  "",
  sprintf("Detected **%d** HGNC genes. Genes with FDR < 0.05: **%d** (%d up, %d down).",
          nrow(tt), summary_dt$fdr_0.05, summary_dt$up_fdr_0.05,
          summary_dt$down_fdr_0.05),
  sprintf("Of these, **%d** also reached |FC| >= 1.5. Genes with unadjusted P < 0.05: %d.",
          summary_dt$fdr_fc_1.5, summary_dt$nominal_p_0.05),
  "",
  "Top five results by raw P value:",
  paste0("- ", top5),
  "",
  "## Methods",
  "",
  "Used the background-corrected, quantile-normalised and log2-transformed matrix provided by the GEO submitters. GPL19615 RefSeq transcripts were mapped to current gene symbols via the frozen HGNC table; the low-expression filter was fixed independently before the primary analysis, and multiple probes per gene were averaged. For each subject the week-4-minus-baseline expression difference was computed, and the mean paired change was tested with a limma robust empirical-Bayes model (trend covariate = mean log2 expression abundance across the 14 primary-analysis samples, fixed 2026-09-10; v1 mistakenly used the row mean of the difference matrix).",
  "",
  "## Interpretation limits",
  "",
  "- The sample is only 7 pairs, so power is limited; leave-one-donor sensitivity is reported alongside.",
  "- Altitude exposure cannot be separated from four weeks of endurance training; results are not a pure hypoxia causal effect.",
  "- 2 260 m is below the 3 800/5 100 m of GSE196728, i.e. a milder adaptation tier.",
  "- GEO provides no AMS scores or individual-level clinical acclimatisation measures.",
  "- Cell composition is captured only by marker-gene proxy scores; with n = 7 no composition covariate adjustment is performed.",
  "",
  "## Key files",
  "",
  "- `differential_expression.tsv`: full gene-level primary results.",
  "- `leave_one_donor_out.tsv`: leave-one-donor stability.",
  "- `probe_mapping_audit.tsv`: platform-to-HGNC mapping audit.",
  "- `cell_composition_proxy_changes.tsv`: cell-composition proxy changes.",
  "- `analysis_registry.json`: frozen contract and interpretation limits."
)
writeLines(summary_lines, file.path(output_dir, "ANALYSIS_SUMMARY.md"),
           useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

artifact_files <- setdiff(
  list.files(output_dir, full.names = TRUE),
  file.path(output_dir, "artifact_manifest.tsv")
)
artifact_manifest <- data.table(
  file = basename(artifact_files),
  bytes = file.info(artifact_files)$size,
  sha256 = vapply(artifact_files, sha256, character(1L))
)
write_tsv(artifact_manifest, file.path(output_dir, "artifact_manifest.tsv"))

message(dataset, " complete: ", ncol(diff), " pairs, ", nrow(tt),
        " genes, ", sum(tt$adj.P.Val < 0.05), " FDR<0.05")
