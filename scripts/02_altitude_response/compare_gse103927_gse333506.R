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
required <- c("gse103927", "gse103927-sensitivity", "gse333506",
              "hgnc", "gene-sets", "output-dir")
missing <- setdiff(required, names(args))
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

analysis_id <- "m2_gse103927_gse333506_agreement_v2_20260910"
g103_path <- normalizePath(args$gse103927, mustWork = TRUE)
g103_sensitivity_path <- normalizePath(args[["gse103927-sensitivity"]],
                                       mustWork = TRUE)
g333_path <- normalizePath(args$gse333506, mustWork = TRUE)
hgnc_path <- normalizePath(args$hgnc, mustWork = TRUE)
gmt_path <- normalizePath(args[["gene-sets"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}
sha256 <- function(path) {
  z <- system2("sha256sum", path, stdout = TRUE)
  sub("  .*", "", z[[1L]])
}

g103 <- fread(g103_path, showProgress = FALSE)
g103_sensitivity <- fread(g103_sensitivity_path, showProgress = FALSE)
g333 <- fread(g333_path, showProgress = FALSE)
if (anyDuplicated(g103$hgnc_id) || anyDuplicated(g333$hgnc_id)) {
  stop("HGNC IDs must be unique within each result")
}

h <- fread(hgnc_path,
           select = c("hgnc_id", "symbol", "status", "locus_group"),
           showProgress = FALSE)
h <- h[status == "Approved"]
common <- merge(
  g103[, .(
    hgnc_id,
    g103_source_symbol = gene_symbol,
    g103_log2fc = logFC,
    g103_t = t,
    g103_p = P.Value,
    g103_fdr = adj.P.Val,
    g103_fdr_sig = adj.P.Val < 0.05,
    g103_fc15_sig = adj.P.Val < 0.05 & abs(logFC) >= log2(1.5)
  )],
  g333[, .(
    hgnc_id,
    g333_source_symbol = gene_symbol,
    g333_log2fc = logFC,
    g333_t = t,
    g333_p = P.Value,
    g333_fdr = adj.P.Val,
    g333_fdr_sig = adj.P.Val < 0.05,
    g333_fc15_sig = adj.P.Val < 0.05 & abs(logFC) >= log2(1.5)
  )],
  by = "hgnc_id", sort = FALSE
)
common <- merge(common, h[, .(hgnc_id, current_symbol = symbol, locus_group)],
                by = "hgnc_id", all.x = TRUE, sort = FALSE)
if (any(is.na(common$current_symbol))) stop("Common HGNC ID absent from registry")
common[, same_direction := sign(g103_log2fc) == sign(g333_log2fc)]
common[, both_fdr := g103_fdr_sig & g333_fdr_sig]
common[, both_fdr_same_direction := both_fdr & same_direction]
common[, both_fc15 := g103_fc15_sig & g333_fc15_sig]
common[, combined_abs_z := abs(g103_t) + abs(g333_t)]
setorder(common, -combined_abs_z)

metric_row <- function(z, stratum) {
  pearson_test <- cor.test(z$g103_log2fc, z$g333_log2fc,
                           method = "pearson")
  spearman_test <- suppressWarnings(cor.test(
    z$g103_log2fc, z$g333_log2fc, method = "spearman", exact = FALSE
  ))
  sign_test <- binom.test(sum(z$same_direction), nrow(z), p = 0.5)
  data.table(
    stratum,
    n_common_genes = nrow(z),
    pearson_log2fc = unname(pearson_test$estimate),
    pearson_p_naive = pearson_test$p.value,
    spearman_log2fc = unname(spearman_test$estimate),
    spearman_p_naive = spearman_test$p.value,
    sign_concordance = mean(z$same_direction),
    sign_concordance_ci_low = sign_test$conf.int[[1L]],
    sign_concordance_ci_high = sign_test$conf.int[[2L]],
    sign_test_p_naive = sign_test$p.value,
    g103_fdr = sum(z$g103_fdr_sig),
    g333_fdr = sum(z$g333_fdr_sig),
    both_fdr = sum(z$both_fdr),
    both_fdr_same_direction = sum(z$both_fdr_same_direction),
    both_fdr_opposite_direction = sum(z$both_fdr & !z$same_direction)
  )
}
agreement <- rbind(
  metric_row(common, "all_common_current_HGNC"),
  metric_row(common[locus_group == "protein-coding gene"], "protein_coding")
)

overlap_test <- function(z, threshold, a_col, b_col) {
  tab <- table(factor(z[[a_col]], levels = c(FALSE, TRUE)),
               factor(z[[b_col]], levels = c(FALSE, TRUE)))
  ft <- fisher.test(tab)
  data.table(
    threshold,
    universe_n = nrow(z),
    g103_positive = sum(z[[a_col]]),
    g333_positive = sum(z[[b_col]]),
    overlap = sum(z[[a_col]] & z[[b_col]]),
    fisher_odds_ratio = unname(ft$estimate),
    fisher_p = ft$p.value
  )
}
overlap <- rbind(
  overlap_test(common, "FDR<0.05", "g103_fdr_sig", "g333_fdr_sig"),
  overlap_test(common, "FDR<0.05_and_absFC>=1.5",
               "g103_fc15_sig", "g333_fc15_sig")
)

conditional <- merge(
  g103_sensitivity[, .(
    hgnc_id,
    g103_log2fc = proxy_conditional_logFC,
    g103_t = proxy_conditional_t,
    g103_p = proxy_conditional_p,
    g103_fdr = proxy_conditional_fdr,
    g103_fdr_sig = proxy_conditional_fdr < 0.05,
    g103_fc15_sig =
      proxy_conditional_fdr < 0.05 &
      abs(proxy_conditional_logFC) >= log2(1.5)
  )],
  g333[, .(
    hgnc_id,
    g333_log2fc = logFC,
    g333_t = t,
    g333_p = P.Value,
    g333_fdr = adj.P.Val,
    g333_fdr_sig = adj.P.Val < 0.05,
    g333_fc15_sig = adj.P.Val < 0.05 & abs(logFC) >= log2(1.5)
  )],
  by = "hgnc_id", sort = FALSE
)
conditional <- merge(
  conditional, h[, .(hgnc_id, current_symbol = symbol, locus_group)],
  by = "hgnc_id", all.x = TRUE, sort = FALSE
)
conditional[, same_direction := sign(g103_log2fc) == sign(g333_log2fc)]
conditional[, both_fdr := g103_fdr_sig & g333_fdr_sig]
conditional[, both_fdr_same_direction := both_fdr & same_direction]
conditional_agreement <- rbind(
  metric_row(conditional, "all_common_current_HGNC"),
  metric_row(conditional[locus_group == "protein-coding gene"],
             "protein_coding")
)

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  pieces <- strsplit(lines, "\t", fixed = TRUE)
  names(pieces) <- vapply(pieces, `[[`, character(1L), 1L)
  lapply(pieces, function(x) unique(x[-c(1L, 2L)]))
}
sets <- read_gmt(gmt_path)
common_by_symbol <- common[!duplicated(current_symbol)]
set_indices <- lapply(sets, function(x) {
  which(common_by_symbol$current_symbol %in% x)
})
set_indices <- set_indices[lengths(set_indices) >= 10L &
                             lengths(set_indices) <= 500L]
if (!length(set_indices)) stop("No eligible pathways in common universe")

run_camera <- function(statistic, prefix) {
  names(statistic) <- common_by_symbol$current_symbol
  z <- as.data.table(cameraPR(statistic, set_indices,
                              inter.gene.cor = 0.01, sort = FALSE),
                     keep.rownames = "pathway")
  setnames(z, c("NGenes", "Direction", "PValue", "FDR"),
           paste0(prefix, c("_n_genes", "_direction", "_p", "_fdr")))
  z[, (paste0(prefix, "_signed_z")) := {
    p <- pmax(get(paste0(prefix, "_p")), .Machine$double.xmin)
    direction <- fifelse(get(paste0(prefix, "_direction")) == "Up", 1, -1)
    direction * qnorm(p / 2, lower.tail = FALSE)
  }]
  z
}
path103 <- run_camera(common_by_symbol$g103_t, "g103")
path333 <- run_camera(common_by_symbol$g333_t, "g333")
pathway <- merge(path103, path333, by = "pathway", sort = FALSE)
pathway[, same_direction := g103_direction == g333_direction]
pathway[, both_fdr := g103_fdr < 0.05 & g333_fdr < 0.05]
pathway[, combined_abs_z := abs(g103_signed_z) + abs(g333_signed_z)]
setorder(pathway, -combined_abs_z)

pathway_pearson <- cor.test(pathway$g103_signed_z, pathway$g333_signed_z,
                            method = "pearson")
pathway_spearman <- suppressWarnings(cor.test(
  pathway$g103_signed_z, pathway$g333_signed_z,
  method = "spearman", exact = FALSE
))
pathway_sign <- binom.test(sum(pathway$same_direction), nrow(pathway), p = 0.5)
pathway_agreement <- data.table(
  n_pathways = nrow(pathway),
  pearson_signed_z = unname(pathway_pearson$estimate),
  pearson_p_naive = pathway_pearson$p.value,
  spearman_signed_z = unname(pathway_spearman$estimate),
  spearman_p_naive = pathway_spearman$p.value,
  sign_concordance = mean(pathway$same_direction),
  sign_concordance_ci_low = pathway_sign$conf.int[[1L]],
  sign_concordance_ci_high = pathway_sign$conf.int[[2L]],
  sign_test_p_naive = pathway_sign$p.value,
  g103_fdr = sum(pathway$g103_fdr < 0.05),
  g333_fdr = sum(pathway$g333_fdr < 0.05),
  both_fdr = sum(pathway$both_fdr),
  both_fdr_same_direction = sum(pathway$both_fdr & pathway$same_direction)
)

source_manifest <- data.table(
  path = c(g103_path, g103_sensitivity_path, g333_path, hgnc_path, gmt_path),
  bytes = file.info(c(g103_path, g103_sensitivity_path, g333_path,
                      hgnc_path, gmt_path))$size,
  sha256 = vapply(c(g103_path, g103_sensitivity_path, g333_path,
                    hgnc_path, gmt_path),
                  sha256, character(1L))
)
registry <- list(
  analysis_id = analysis_id,
  status = "complete",
  completed_at = "2026-07-28",
  estimand_gse103927 = "PBMC ALT16 at 5260 m minus sea-level baseline",
  estimand_gse333506 =
    "whole-blood Han athlete week 4 at 2260 m minus pre-ascent baseline",
  direction = "positive values mean higher expression after altitude exposure",
  gene_alignment = "exact frozen current HGNC ID intersection",
  pooling_rule = "directional comparison only; no pooled effect",
  pathway_method =
    "cameraPR on study-specific moderated t statistics over identical common-gene universe",
  composition_sensitivity =
    "GSE103927 paired differences regressed on four uncentered PBMC marker-score changes; intercept estimates zero-proxy-change response",
  inference_boundary =
    "reported correlation and sign-test P values are naive because genes and pathways are dependent",
  key_heterogeneity =
    "PBMC versus whole blood; 5260 m expedition versus 2260 m altitude training; day 16 versus week 4"
)

write_tsv(common, file.path(output_dir, "gene_direction_comparison.tsv"))
write_tsv(agreement, file.path(output_dir, "gene_direction_agreement.tsv"))
write_tsv(overlap, file.path(output_dir, "significant_gene_overlap.tsv"))
write_tsv(conditional,
          file.path(output_dir,
                    "composition_proxy_gene_direction_comparison.tsv"))
write_tsv(conditional_agreement,
          file.path(output_dir,
                    "composition_proxy_direction_agreement.tsv"))
write_tsv(common[same_direction == TRUE][1:min(.N, 100L)],
          file.path(output_dir, "top_concordant_genes.tsv"))
write_tsv(common[same_direction == FALSE][1:min(.N, 100L)],
          file.path(output_dir, "top_discordant_genes.tsv"))
write_tsv(pathway, file.path(output_dir, "pathway_direction_comparison.tsv"))
write_tsv(pathway_agreement,
          file.path(output_dir, "pathway_direction_agreement.tsv"))
write_tsv(source_manifest, file.path(output_dir, "source_manifest.tsv"))
write_json(registry, file.path(output_dir, "analysis_registry.json"),
           pretty = TRUE, auto_unbox = TRUE)

scatter <- ggplot(common, aes(g333_log2fc, g103_log2fc,
                              color = same_direction)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_point(alpha = 0.35, size = 0.8) +
  scale_color_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B")) +
  labs(x = "GSE333506 log2FC (week 4 - baseline)",
       y = "GSE103927 log2FC (ALT16 - sea level)",
       title = "Altitude-response directions in common HGNC genes",
       color = "same direction") +
  theme_bw(base_size = 11)
ggsave(file.path(output_dir, "gene_effect_scatter.png"), scatter,
       width = 7, height = 6, dpi = 180)

path_scatter <- ggplot(pathway, aes(g333_signed_z, g103_signed_z,
                                    color = same_direction)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_point(alpha = 0.55, size = 1.2) +
  scale_color_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B")) +
  labs(x = "GSE333506 pathway signed Z",
       y = "GSE103927 pathway signed Z",
       title = "Pathway directions on the common gene universe",
       color = "same direction") +
  theme_bw(base_size = 11)
ggsave(file.path(output_dir, "pathway_effect_scatter.png"), path_scatter,
       width = 7, height = 6, dpi = 180)

a <- agreement[stratum == "all_common_current_HGNC"]
ac <- conditional_agreement[stratum == "all_common_current_HGNC"]
pway <- pathway_agreement[1L]
summary_lines <- c(
  "# GSE103927 vs GSE333506 directional-agreement summary",
  "",
  "## Comparisons",
  "",
  "- GSE103927: PBMC of 21 healthy lowland subjects, day 16 at 5 260 m minus sea-level baseline.",
  "- GSE333506: whole blood of 7 Han Chinese male athletes, week 4 at about 2 260 m minus pre-altitude baseline.",
  "- For both, positive logFC means higher expression after altitude exposure; the intersection uses frozen current HGNC IDs.",
  "",
  "## Gene-level results",
  "",
  sprintf("Common genes: **%d**; log2FC Pearson r = %.3f, Spearman rho = %.3f.",
          a$n_common_genes, a$pearson_log2fc, a$spearman_log2fc),
  sprintf("Sign concordance: **%.1f%%** (95%% binomial interval %.1f%%-%.1f%%).",
          100 * a$sign_concordance, 100 * a$sign_concordance_ci_low,
          100 * a$sign_concordance_ci_high),
  sprintf("Among common genes, %d are FDR-significant in GSE103927 and %d in GSE333506; %d are significant in both (%d same direction, %d opposite).",
          a$g103_fdr, a$g333_fdr, a$both_fdr,
          a$both_fdr_same_direction, a$both_fdr_opposite_direction),
  sprintf("After composition-proxy conditioning, GSE103927 has %d FDR-significant genes; versus GSE333506: Pearson r = %.3f, Spearman rho = %.3f, sign concordance %.1f%%.",
          ac$g103_fdr, ac$pearson_log2fc, ac$spearman_log2fc,
          100 * ac$sign_concordance),
  "",
  "## Pathway-level results",
  "",
  sprintf("%d pathways tested on the same common gene universe; signed-Z Pearson r = %.3f, Spearman rho = %.3f, direction concordance %.1f%%.",
          pway$n_pathways, pway$pearson_signed_z, pway$spearman_signed_z,
          100 * pway$sign_concordance),
  sprintf("%d pathways are FDR-significant in GSE103927 and %d in GSE333506; %d in both (%d same direction).",
          pway$g103_fdr, pway$g333_fdr, pway$both_fdr,
          pway$both_fdr_same_direction),
  "",
  "## Interpretation limits",
  "",
  "- This is a directional-replication analysis; no combined effect sizes are computed.",
  "- PBMC vs whole blood, a 5 260 m expedition vs 2 260 m training, and 16 days vs 4 weeks are all prespecified sources of heterogeneity.",
  "- P values of the correlation and sign tests are descriptive only, because genes and pathways are not independent.",
  "- Whether the results constitute biological replication should be judged from effect correlations, sign concordance and pathway results together, not from the intersection of significant genes alone."
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

message("Comparison complete: ", nrow(common), " common genes, ",
        sprintf("r=%.3f, sign=%.1f%%", a$pearson_log2fc,
                100 * a$sign_concordance))
