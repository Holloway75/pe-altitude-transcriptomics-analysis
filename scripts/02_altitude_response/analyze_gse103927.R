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

dataset <- "GSE103927"
analysis_id <- "m2_gse103927_alt16_paired_v2_20260910"
geo_dir <- normalizePath(args[["geo-dir"]], mustWork = TRUE)
hgnc_path <- normalizePath(args[["hgnc"]], mustWork = TRUE)
frozen_dir <- normalizePath(args[["frozen-dir"]], mustWork = FALSE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
dir.create(frozen_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(geo_dir, "GSE103927_series_matrix.txt.gz")
platform_path <- file.path(geo_dir, "GPL6244.annot.gz")
stopifnot(file.exists(matrix_path), file.exists(platform_path))

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
    if (length(z) != 1L) stop("Expected one metadata line: ", prefix)
    quoted_values(z)
  }
  gsm <- get_one("!Sample_geo_accession")
  meta <- data.table(
    sample_id = gsm,
    title = get_one("!Sample_title"),
    source_name = get_one("!Sample_source_name_ch1")
  )
  characteristics <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
  for (line in characteristics) {
    values <- quoted_values(line)
    if (length(values) != length(gsm)) stop("Characteristic length mismatch")
    key <- sub(":.*$", "", values[[1L]])
    value <- trimws(sub("^[^:]+:", "", values))
    if (key %in% names(meta)) stop("Duplicate characteristic key: ", key)
    meta[, (key) := value]
  }
  required_fields <- c("time point", "subject", "cell type", "age", "Sex")
  if (length(setdiff(required_fields, names(meta)))) {
    stop("Missing sample characteristics: ",
         paste(setdiff(required_fields, names(meta)), collapse = ", "))
  }
  setnames(meta, c("time point", "subject", "cell type", "Sex"),
           c("timepoint", "subject_id", "cell_type", "sex"))
  meta
}

read_expression <- function(path) {
  x <- fread(path, skip = "!series_matrix_table_begin", header = TRUE,
             showProgress = FALSE)
  if (ncol(x) != 113L) stop("Expected ID_REF plus 112 expression columns")
  feature_id <- as.character(x[[1L]])
  x[[1L]] <- NULL
  expr <- as.matrix(x)
  storage.mode(expr) <- "double"
  rownames(expr) <- feature_id
  if (any(!is.finite(expr))) stop("Non-finite expression values")
  expr
}

build_symbol_map <- function(h) {
  exact <- unique(h[, .(source_symbol = symbol, current_symbol = symbol,
                         hgnc_id, mapping_type = "approved_exact")])
  expand_pipe <- function(column, type) {
    z <- h[!is.na(get(column)) & get(column) != "",
           .(source_symbol = unlist(strsplit(get(column), "\\|", perl = TRUE))),
           by = .(current_symbol = symbol, hgnc_id)]
    z[, mapping_type := type]
    z
  }
  alias <- rbind(expand_pipe("prev_symbol", "previous_symbol"),
                 expand_pipe("alias_symbol", "alias_symbol"))
  alias <- alias[!source_symbol %in% exact$source_symbol]
  unique_alias <- alias[, .(n_current = uniqueN(current_symbol)),
                        by = source_symbol][n_current == 1L, source_symbol]
  rbind(exact, unique(alias[source_symbol %in% unique_alias],
                      by = "source_symbol"), fill = TRUE)
}

read_platform_mapping <- function(path, hgnc_path) {
  p <- fread(path, skip = "!platform_table_begin", header = TRUE,
             select = c("ID", "Gene title", "Gene symbol", "Gene ID"),
             showProgress = FALSE)
  p <- p[ID != "!platform_table_end"]
  setnames(p, c("ID", "Gene title", "Gene symbol", "Gene ID"),
           c("feature_id", "gene_title_raw", "gene_symbol_raw",
             "entrez_id_raw"))
  p[, feature_id := as.character(feature_id)]

  h <- fread(hgnc_path,
             select = c("hgnc_id", "symbol", "name", "status", "locus_group",
                        "entrez_id", "alias_symbol", "prev_symbol"),
             showProgress = FALSE)
  h <- h[status == "Approved"]
  h[, entrez_id := as.character(entrez_id)]
  entrez_counts <- h[!is.na(entrez_id) & entrez_id != "",
                     .(n_current = uniqueN(symbol)), by = entrez_id]
  valid_entrez <- entrez_counts[n_current == 1L, entrez_id]
  entrez_map <- unique(
    h[entrez_id %in% valid_entrez,
      .(entrez_id, hgnc_id, current_symbol = symbol,
        gene_description = name, locus_group)],
    by = "entrez_id"
  )
  symbol_map <- build_symbol_map(h)

  resolve_entrez <- function(value) {
    ids <- trimws(strsplit(ifelse(is.na(value), "", value),
                           "///", fixed = TRUE)[[1L]])
    ids <- ids[ids != ""]
    rows <- entrez_map[match(ids, entrez_id), ]
    rows <- unique(rows[!is.na(current_symbol)], by = "current_symbol")
    if (nrow(rows) == 1L) {
      return(list(n = 1L, hgnc_id = rows$hgnc_id,
                  current_symbol = rows$current_symbol,
                  gene_description = rows$gene_description,
                  locus_group = rows$locus_group))
    }
    list(n = nrow(rows), hgnc_id = NA_character_,
         current_symbol = NA_character_, gene_description = NA_character_,
         locus_group = NA_character_)
  }
  resolved <- lapply(p$entrez_id_raw, resolve_entrez)
  p[, n_entrez_current := vapply(resolved, `[[`, integer(1L), "n")]
  p[, hgnc_id := vapply(resolved, `[[`, character(1L), "hgnc_id")]
  p[, gene_symbol := vapply(resolved, `[[`, character(1L), "current_symbol")]
  p[, gene_description :=
      vapply(resolved, `[[`, character(1L), "gene_description")]
  p[, locus_group := vapply(resolved, `[[`, character(1L), "locus_group")]
  p[, mapping_method := fifelse(n_entrez_current == 1L,
                                "unique_current_HGNC_entrez",
                         fifelse(n_entrez_current > 1L,
                                 "ambiguous_multi_gene_entrez",
                                 "unmapped_entrez"))]

  fallback_i <- which(p$mapping_method == "unmapped_entrez" &
                        !is.na(p$gene_symbol_raw) & p$gene_symbol_raw != "")
  for (i in fallback_i) {
    symbols <- trimws(strsplit(p$gene_symbol_raw[[i]], "///",
                               fixed = TRUE)[[1L]])
    rows <- unique(symbol_map[match(symbols, source_symbol), ],
                   by = "current_symbol")
    rows <- rows[!is.na(current_symbol)]
    if (nrow(rows) == 1L) {
      hrow <- h[symbol == rows$current_symbol[[1L]]][1L]
      p[i, `:=`(
        hgnc_id = hrow$hgnc_id,
        gene_symbol = hrow$symbol,
        gene_description = hrow$name,
        locus_group = hrow$locus_group,
        mapping_method = paste0("unique_current_HGNC_",
                                rows$mapping_type[[1L]])
      )]
    } else if (nrow(rows) > 1L) {
      p[i, mapping_method := "ambiguous_multi_gene_symbol"]
    }
  }
  p
}

meta <- read_series_metadata(matrix_path)
expr_all <- read_expression(matrix_path)
if (!identical(colnames(expr_all), meta$sample_id)) {
  stop("Series-matrix sample order does not match metadata")
}

meta[, include_primary := timepoint %in% c("baseline", "day 16")]
meta[, exclusion_reason := fifelse(
  include_primary, NA_character_, paste0(timepoint, " excluded from ALT16 contrast")
)]
primary_meta <- meta[include_primary == TRUE]
pair_check <- primary_meta[, .(
  n = .N,
  n_baseline = sum(timepoint == "baseline"),
  n_alt16 = sum(timepoint == "day 16")
), by = subject_id]
if (nrow(pair_check) != 21L ||
    any(pair_check$n != 2L | pair_check$n_baseline != 1L |
        pair_check$n_alt16 != 1L)) {
  stop("Primary pairing contract failed")
}

baseline <- primary_meta[timepoint == "baseline"][order(as.integer(subject_id))]
alt16 <- primary_meta[timepoint == "day 16"][order(as.integer(subject_id))]
if (!identical(baseline$subject_id, alt16$subject_id)) stop("Pair order mismatch")

platform <- read_platform_mapping(platform_path, hgnc_path)
ann <- platform[match(rownames(expr_all), feature_id)]
if (any(is.na(ann$feature_id))) stop("Matrix probes absent from GPL annotation")
keep_probe <- !is.na(ann$gene_symbol) & ann$gene_symbol != "" &
  grepl("^unique_current_HGNC", ann$mapping_method)
if (sum(keep_probe) < 10000L) stop("Too few uniquely mapped probes")

expr_probe <- expr_all[keep_probe, primary_meta$sample_id, drop = FALSE]
ann_probe <- ann[keep_probe]
expr_gene <- avereps(expr_probe, ID = ann_probe$gene_symbol)
gene_ann <- ann_probe[, .(
  hgnc_id = hgnc_id[[1L]],
  gene_description = gene_description[[1L]],
  locus_group = locus_group[[1L]],
  n_retained_probes = .N,
  retained_probe_ids = paste(feature_id, collapse = "|"),
  mapping_methods = paste(sort(unique(mapping_method)), collapse = "|")
), by = gene_symbol]
gene_ann <- gene_ann[match(rownames(expr_gene), gene_symbol)]

diff <- expr_gene[, alt16$sample_id, drop = FALSE] -
  expr_gene[, baseline$sample_id, drop = FALSE]
colnames(diff) <- baseline$subject_id
# 2026-09-10 fix: the trend covariate must be gene-expression abundance
# (mean log2 intensity across all 42 primary-analysis samples), not the row
# mean of the difference matrix (under an intercept-only design the latter
# equals logFC itself, i.e. building the prior from the effect; see
# docs/report/PRE_SUBMISSION_CODE_AUDIT_20260910.md, finding A-1).
stopifnot(identical(rownames(diff), rownames(expr_gene)))
abundance <- unname(rowMeans(expr_gene))
design <- matrix(1, nrow = ncol(diff), ncol = 1,
                 dimnames = list(colnames(diff), "alt16_minus_sea_level"))
fit <- eBayes(lmFit(diff, design), trend = abundance, robust = TRUE)
tt <- as.data.table(topTable(fit, coef = 1, number = Inf, sort.by = "none"),
                    keep.rownames = "gene_symbol")
tt[, AveExpr := NULL]
tt <- merge(gene_ann, tt, by = "gene_symbol", all.y = TRUE, sort = FALSE)
tt[, mean_expression :=
     rowMeans(expr_gene)[match(gene_symbol, rownames(expr_gene))]]
tt[, direction := fifelse(logFC > 0, "up", "down")]
tt[, significant_fdr_0.05 := adj.P.Val < 0.05]
tt[, significant_fdr_fc_1.5 :=
     adj.P.Val < 0.05 & abs(logFC) >= log2(1.5)]
setcolorder(tt, c("gene_symbol", "hgnc_id", "gene_description", "locus_group",
                  "n_retained_probes", "retained_probe_ids", "mapping_methods",
                  "logFC", "mean_expression", "t", "P.Value", "adj.P.Val", "B",
                  "direction", "significant_fdr_0.05",
                  "significant_fdr_fc_1.5"))
setorder(tt, P.Value)

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
  NK_cell = c("NKG7", "GNLY", "KLRD1", "PRF1")
)
scores <- rbindlist(lapply(names(marker_sets), function(s) {
  k <- intersect(marker_sets[[s]], rownames(expr_gene))
  if (length(k) < 2L) return(NULL)
  data.table(sample_id = colnames(expr_gene), cell_proxy = s,
             score = colMeans(expr_gene[k, , drop = FALSE]),
             marker_count = length(k))
}), use.names = TRUE)
scores <- merge(scores, primary_meta[, .(sample_id, subject_id, timepoint)],
                by = "sample_id")
score_wide <- dcast(scores, subject_id + timepoint ~ cell_proxy,
                    value.var = "score")
score_base <- score_wide[timepoint == "baseline"][order(as.integer(subject_id))]
score_alt16 <- score_wide[timepoint == "day 16"][order(as.integer(subject_id))]
proxy_names <- setdiff(names(score_wide), c("subject_id", "timepoint"))
proxy_change <- as.data.table(
  as.matrix(score_alt16[, ..proxy_names] - score_base[, ..proxy_names])
)
proxy_change[, subject_id := score_base$subject_id]
setcolorder(proxy_change, c("subject_id", proxy_names))

proxy_model <- proxy_change[match(colnames(diff), subject_id)]
proxy_design <- cbind(
  no_proxy_change_intercept = 1,
  as.matrix(proxy_model[, ..proxy_names])
)
colnames(proxy_design)[-1L] <- paste0("delta_", proxy_names)
if (qr(proxy_design)$rank != ncol(proxy_design)) {
  stop("Composition-proxy sensitivity design is not full rank")
}
fit_proxy <- eBayes(lmFit(diff, proxy_design), trend = abundance, robust = TRUE)
tt_proxy <- as.data.table(
  topTable(fit_proxy, coef = "no_proxy_change_intercept",
           number = Inf, sort.by = "none"),
  keep.rownames = "gene_symbol"
)
tt_proxy[, AveExpr := NULL]
setnames(tt_proxy,
         c("logFC", "t", "P.Value", "adj.P.Val", "B"),
         c("proxy_conditional_logFC", "proxy_conditional_t",
           "proxy_conditional_p", "proxy_conditional_fdr",
           "proxy_conditional_B"))
proxy_sensitivity <- merge(
  tt[, .(gene_symbol, hgnc_id, primary_logFC = logFC, primary_t = t,
         primary_p = P.Value, primary_fdr = adj.P.Val)],
  tt_proxy, by = "gene_symbol", sort = FALSE
)
proxy_sensitivity[, same_direction :=
                    sign(primary_logFC) == sign(proxy_conditional_logFC)]
proxy_sensitivity[, proxy_conditional_fdr_sig :=
                    proxy_conditional_fdr < 0.05]
setorder(proxy_sensitivity, proxy_conditional_p)

proxy_summary <- rbindlist(lapply(proxy_names, function(n) {
  x <- proxy_change[[n]]
  w <- suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE))
  data.table(
    cell_proxy = n,
    n_pairs = length(x),
    mean_change = mean(x),
    median_change = median(x),
    wilcoxon_p = w$p.value
  )
}))
proxy_summary[, wilcoxon_fdr := p.adjust(wilcoxon_p, method = "BH")]

sample_out <- meta[, .(
  dataset_id = dataset,
  sample_id,
  expression_column = sample_id,
  donor_id = subject_id,
  timepoint,
  include_primary,
  exclusion_reason,
  age = as.integer(age),
  sex,
  tissue = "peripheral blood mononuclear cells",
  assay = "GPL6244 Affymetrix Human Gene 1.0 ST array; GEO RMA log2",
  altitude_m = fifelse(timepoint %in% c("baseline", "day 7 post-descent",
                                        "day 21 post-descent"),
                       fifelse(timepoint == "baseline", 0L, 1525L), 5260L),
  visit = timepoint,
  pair_id = subject_id,
  source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103927",
  geo_title = title,
  source_name,
  submitted_cell_type = cell_type
)]
pair_out <- data.table(
  donor_id = baseline$subject_id,
  baseline_sample = baseline$sample_id,
  baseline_title = baseline$title,
  alt16_sample = alt16$sample_id,
  alt16_title = alt16$title
)
mapping_audit <- platform[, .(
  n_probes = .N,
  n_in_expression_matrix = sum(feature_id %in% rownames(expr_all)),
  n_retained = sum(feature_id %in% ann_probe$feature_id)
), by = mapping_method]
setorder(mapping_audit, mapping_method)

source_paths <- c(matrix_path, platform_path, hgnc_path)
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
    "abundance across the 42 primary samples; see",
    "docs/report/PRE_SUBMISSION_CODE_AUDIT_20260910.md (finding A-1).",
    "Donor difference matrices, logFC and ordinary one-sample statistics are",
    "unchanged relative to v1; moderated t/P/FDR shift modestly."
  ),
  population = "21 healthy lowlanders (9 women, 12 men)",
  n_samples_total = nrow(meta),
  n_samples_primary = nrow(primary_meta),
  n_donors_primary = nrow(pair_check),
  primary_contrast = "day 16 at 5260 m minus sea-level baseline",
  direction = "positive logFC means higher expression at ALT16",
  tissue = "PBMC",
  value_type = "GEO-submitted RMA-normalized log2 intensity",
  feature_mapping =
    "GPL6244 Entrez IDs, with unique-symbol fallback, to approved frozen HGNC",
  independent_filter =
    "unique current-HGNC mapping only; no intensity threshold because detection calls are unavailable",
  probe_collapse = "mean RMA intensity across retained probes per HGNC symbol",
  primary_model = paste(
    "one-sample empirical-Bayes limma (robust, trend on mean log2 expression",
    "abundance across the 42 primary samples) on 21 within-donor",
    "ALT16-baseline differences"
  ),
  multiple_testing = "Benjamini-Hochberg",
  significance_threshold = "FDR < 0.05",
  composition_rule = "PBMC marker-score diagnostics only; no covariate adjustment",
  confirmatory_boundary =
    "PBMC rather than whole blood; expedition exposure at 5260 m; not directly poolable with GSE333506"
)

write_tsv(sample_out, file.path(frozen_dir, "sample_metadata.tsv"))
write_tsv(pair_out, file.path(frozen_dir, "donor_altitude_mapping.tsv"))
write_tsv(mapping_audit, file.path(frozen_dir, "probe_mapping_audit.tsv"))
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
write_tsv(proxy_summary,
          file.path(output_dir, "cell_composition_proxy_summary.tsv"))
write_tsv(proxy_sensitivity,
          file.path(output_dir, "composition_proxy_sensitivity.tsv"))
write_tsv(sample_out, file.path(output_dir, "sample_metadata.tsv"))
write_tsv(pair_out, file.path(output_dir, "donor_altitude_mapping.tsv"))
write_tsv(mapping_audit, file.path(output_dir, "probe_mapping_audit.tsv"))
write_tsv(source_manifest, file.path(output_dir, "source_manifest.tsv"))
write_json(c(freeze, list(analysis_status = "complete")),
           file.path(output_dir, "analysis_registry.json"),
           pretty = TRUE, auto_unbox = TRUE)

summary_dt <- data.table(
  dataset_id = dataset,
  contrast = "alt16_5260m_minus_sea_level_baseline",
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
    cell_proxy_changes = proxy_change,
    composition_proxy_design = proxy_design,
    composition_proxy_fit = fit_proxy
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
  cor(expr_gene[, baseline$sample_id[i]], expr_gene[, alt16$sample_id[i]])
}, numeric(1L))
pair_cor_dt <- rbind(
  data.table(sample_id = baseline$sample_id,
             within_pair_correlation = pair_cor),
  data.table(sample_id = alt16$sample_id,
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
       title = "GSE103927 primary-analysis sample PCA") +
  theme_bw(base_size = 11)
ggsave(file.path(output_dir, "pca.png"), p, width = 7, height = 5, dpi = 160)

box_dt <- as.data.table(expr_gene, keep.rownames = "gene_symbol")
box_dt <- melt(box_dt, id.vars = "gene_symbol",
               variable.name = "sample_id", value.name = "expression")
box_dt <- merge(box_dt,
                primary_meta[, .(sample_id, subject_id, timepoint)],
                by = "sample_id")
pb <- ggplot(box_dt, aes(sample_id, expression, fill = timepoint)) +
  geom_boxplot(outlier.size = 0.15) +
  labs(x = NULL, y = "GEO RMA log2 intensity",
       title = "GSE103927 primary-analysis expression distribution") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
ggsave(file.path(output_dir, "normalized_expression_boxplot.png"),
       pb, width = 12, height = 5, dpi = 160)

volcano <- copy(tt)
volcano[, neglog10p := -log10(pmax(P.Value, .Machine$double.xmin))]
pv <- ggplot(volcano, aes(logFC, neglog10p,
                          color = significant_fdr_0.05)) +
  geom_point(alpha = 0.55, size = 1) +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#777777")) +
  labs(x = "Paired log2 FC (ALT16 - sea-level baseline)", y = "-log10(P)",
       title = "GSE103927 paired differential expression") +
  theme_bw(base_size = 11) + guides(color = "none")
ggsave(file.path(output_dir, "volcano.png"), pv, width = 7, height = 5,
       dpi = 160)

top5 <- tt[1:min(5L, .N),
           paste0(gene_symbol, " (logFC=", sprintf("%.3f", logFC),
                  ", FDR=", format(adj.P.Val, digits = 3), ")")]
summary_lines <- c(
  "# GSE103927 ALT16 paired differential-expression summary",
  "",
  "## Primary estimand",
  "",
  "Included 21 healthy lowland subjects; the primary comparison is day 16 of acclimatisation at 5 260 m minus the sea-level baseline within the same subject. Other time points do not enter the primary model.",
  "",
  "## Results",
  "",
  sprintf("Detected **%d** current-HGNC genes. Genes with FDR < 0.05: **%d** (%d up, %d down).",
          nrow(tt), summary_dt$fdr_0.05, summary_dt$up_fdr_0.05,
          summary_dt$down_fdr_0.05),
  sprintf("Of these, **%d** also reached |FC| >= 1.5; genes with unadjusted P < 0.05: %d.",
          summary_dt$fdr_fc_1.5, summary_dt$nominal_p_0.05),
  "",
  "Top five results by raw P value:",
  paste0("- ", top5),
  "",
  "## Methods",
  "",
  "Used the RMA-normalised log2 matrix provided by the GEO submitters. GPL6244 probes were mapped via Entrez ID to the frozen current-HGNC records where possible, keeping only unique mappings and averaging multiple probes per gene. For each subject the ALT16-minus-sea-level-baseline expression difference was computed, and the mean paired change was tested with a limma robust empirical-Bayes model (trend covariate = mean log2 expression abundance across the 42 primary-analysis samples, fixed 2026-09-10; v1 mistakenly used the row mean of the difference matrix).",
  "",
  "## Interpretation limits",
  "",
  "- GSE103927 profiles PBMC whereas GSE333506 profiles whole blood; cross-study agreement is replication across blood compartments.",
  "- The exposure is a 5 260 m expedition setting; altitude cannot be separated from travel and other accompanying environmental factors.",
  "- The processed matrix carries no detection P values, so filtering used unique-HGNC mapping only, with no arbitrary intensity threshold.",
  "- Cell composition was assessed only diagnostically via PBMC marker-gene proxy scores and does not enter the primary model.",
  sprintf("- The composition-proxy-conditional sensitivity model has %d genes at FDR < 0.05; its intercept is the extrapolated effect when all four proxies are unchanged and does not replace the primary model.",
          sum(proxy_sensitivity$proxy_conditional_fdr_sig)),
  "",
  "## Key files",
  "",
  "- `differential_expression.tsv`: full gene-level primary results.",
  "- `leave_one_donor_out.tsv`: leave-one-donor stability.",
  "- `probe_mapping_audit.tsv`: GPL6244-to-current-HGNC mapping audit.",
  "- `cell_composition_proxy_changes.tsv`: PBMC composition proxy changes.",
  "- `composition_proxy_sensitivity.tsv`: composition-proxy-conditional sensitivity results.",
  "- `analysis_registry.json`: frozen estimand and interpretation limits."
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
