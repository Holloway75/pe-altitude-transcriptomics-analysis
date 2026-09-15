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

read_gmt <- function(path) {
  pieces <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  names(pieces) <- vapply(pieces, `[[`, character(1L), 1L)
  lapply(pieces, function(x) unique(x[-c(1L, 2L)]))
}

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)

args <- parse_args(commandArgs(trailingOnly = TRUE))
config_path <- normalizePath(args[["config"]], mustWork = TRUE)
locked_dir <- normalizePath(args[["locked-dir"]], mustWork = TRUE)
result_dir <- normalizePath(args[["result-dir"]], mustWork = TRUE)
module2_release <- normalizePath(args[["module2-release"]], mustWork = TRUE)
hallmark_path <- normalizePath(args[["hallmark-gmt"]], mustWork = TRUE)
reactome_path <- normalizePath(args[["reactome-gmt"]], mustWork = TRUE)

config <- fromJSON(config_path, simplifyVector = TRUE)
universe_path <- file.path(locked_dir, "gene_universe_and_mapping.tsv")
registry_path <- file.path(locked_dir, "module3_v5_analysis_registry.json")
if (!file.exists(universe_path) || !file.exists(registry_path)) {
  stop("Overlap phase must complete before compact signatures")
}
signature_output <- file.path(result_dir, "hypoxia_altitude_signatures.tsv")
core_output <- file.path(result_dir, "signature_core_genes.tsv")
signature_registry_output <- file.path(locked_dir, "compact_signature_registry.tsv")
if (any(file.exists(c(signature_output, core_output, signature_registry_output)))) {
  stop("Refusing to overwrite existing compact-signature output")
}

universe <- fread(universe_path)
if (nrow(universe) != 8101L || anyDuplicated(universe$hgnc_id)) {
  stop("Expected unique 8,101-gene universe")
}

formal_path <- file.path(
  module2_release, "donor_contrasts/ALT_GSE103927_ALT16_TOTAL.tsv.gz"
)
context_path <- file.path(
  module2_release, "donor_contrasts/ALT_GSE333506_HAN_WEEK4_TOTAL.tsv.gz"
)
formal <- fread(formal_path, header = TRUE)
context <- fread(context_path, header = TRUE)
formal <- formal[match(universe$hgnc_id, hgnc_id)]
if (anyNA(formal$hgnc_id)) stop("Formal donor matrix does not cover universe")

donor_columns <- setdiff(names(formal), c("hgnc_id", "hgnc_symbol"))
formal_matrix <- as.matrix(formal[, ..donor_columns])
storage.mode(formal_matrix) <- "double"
rownames(formal_matrix) <- universe$hgnc_symbol
if (ncol(formal_matrix) != 21L || any(!is.finite(formal_matrix))) {
  stop("Invalid formal paired-difference matrix")
}

all_sets <- c(read_gmt(hallmark_path), read_gmt(reactome_path))
requested <- unname(config$compact_signatures)
missing <- setdiff(requested, names(all_sets))
if (length(missing)) stop("Missing frozen signatures: ", paste(missing, collapse = ", "))
sets <- all_sets[requested]
indices <- lapply(sets, function(genes) which(rownames(formal_matrix) %in% genes))
if (any(lengths(indices) < 5L)) stop("One or more signatures has fewer than 5 universe genes")

signature_registry <- rbindlist(lapply(seq_along(sets), function(i) {
  source_path <- if (grepl("^HALLMARK_", names(sets)[i])) hallmark_path else reactome_path
  source_genes <- sets[[i]]
  retained_genes <- sort(intersect(source_genes, rownames(formal_matrix)))
  data.table(
    analysis_id = config$analysis_id,
    signature_id = names(sets)[i],
    source = if (grepl("^HALLMARK_", names(sets)[i])) "MSigDB_2024.1_Hallmark" else "MSigDB_2024.1_Reactome",
    source_path = source_path,
    source_sha256 = sha256(source_path),
    source_n = length(source_genes),
    universe_n = length(retained_genes),
    universe_genes = paste(retained_genes, collapse = "|")
  )
}))
fwrite(signature_registry, signature_registry_output, sep = "\t", quote = FALSE)

design <- matrix(1, nrow = ncol(formal_matrix), ncol = 1L,
                 dimnames = list(colnames(formal_matrix), "mean_difference"))
camera_result <- as.data.table(camera(
  formal_matrix,
  index = indices,
  design = design,
  contrast = 1,
  inter.gene.cor = 0.01,
  use.ranks = FALSE,
  sort = FALSE,
  directional = TRUE
), keep.rownames = "signature_id")
setnames(camera_result, c("NGenes", "Direction", "PValue", "FDR"),
         c("formal_n_genes", "formal_direction", "formal_p", "formal_bh_fdr"))

# 2026-09-10: take the module-2 frozen moderated t directly (consistent with
# the primary model). An earlier version re-fitted in place without trend,
# inconsistent with the module-2 primary model (audit finding Minor M3-6).
formal_t <- universe$alt_module2_t
names(formal_t) <- universe$hgnc_symbol
if (anyNA(formal_t) || any(!is.finite(formal_t))) {
  stop("Universe is missing finite module-2 moderated t values")
}
core_rows <- list()
core_text <- character(length(sets))
for (i in seq_along(sets)) {
  signature <- names(sets)[i]
  member_t <- formal_t[indices[[i]]]
  direction <- camera_result[signature_id == signature, formal_direction]
  aligned <- if (direction == "Up") sort(member_t, decreasing = TRUE) else sort(member_t)
  aligned <- head(aligned, 10L)
  core_text[i] <- paste(names(aligned), collapse = "|")
  core_rows[[i]] <- data.table(
    analysis_id = config$analysis_id,
    signature_id = signature,
    rank = seq_along(aligned),
    hgnc_symbol = names(aligned),
    moderated_t = unname(aligned),
    aligned_with_signature_direction = TRUE,
    role = "descriptive_top_aligned_members_not_leading_edge_inference"
  )
}
camera_result[, core_contributing_genes := core_text[match(signature_id, names(sets))]]

significant <- camera_result[formal_bh_fdr < 0.05, signature_id]
camera_result[, `:=`(
  context_status = "NOT_TRIGGERED_PRIMARY_SIGNATURE_FDR_GE_0.05",
  context_n_genes = NA_integer_,
  context_direction = NA_character_,
  context_p_descriptive = NA_real_,
  context_bh_fdr_descriptive = NA_real_
)]
if (length(significant)) {
  context_symbols <- context$hgnc_symbol
  context_donor_columns <- setdiff(names(context), c("hgnc_id", "hgnc_symbol"))
  context_matrix <- as.matrix(context[, ..context_donor_columns])
  storage.mode(context_matrix) <- "double"
  rownames(context_matrix) <- context_symbols
  context_sets <- sets[significant]
  context_indices <- lapply(context_sets, function(genes) which(context_symbols %in% genes))
  eligible <- lengths(context_indices) >= 5L
  if (any(eligible)) {
    context_design <- matrix(1, nrow = ncol(context_matrix), ncol = 1L)
    context_camera <- as.data.table(camera(
      context_matrix,
      index = context_indices[eligible],
      design = context_design,
      contrast = 1,
      inter.gene.cor = 0.01,
      use.ranks = FALSE,
      sort = FALSE,
      directional = TRUE
    ), keep.rownames = "signature_id")
    setnames(context_camera, c("NGenes", "Direction", "PValue", "FDR"),
             c("context_n_genes", "context_direction", "context_p_descriptive",
               "context_bh_fdr_descriptive"))
    camera_result[context_camera, on = "signature_id", `:=`(
      context_status = "TRIGGERED_DIRECTIONAL_CONTEXT_NOT_REPLICATION",
      context_n_genes = i.context_n_genes,
      context_direction = i.context_direction,
      context_p_descriptive = i.context_p_descriptive,
      context_bh_fdr_descriptive = i.context_bh_fdr_descriptive
    )]
  }
}

camera_result[, `:=`(
  analysis_id = config$analysis_id,
  inference_role = "altitude_axis_biological_positive_control_not_PE_overlap_test",
  formal_axis = "GSE103927_ALT16_minus_baseline_paired_difference",
  camera_inter_gene_correlation = 0.01,
  multiplicity_family_n = length(sets)
)]
setcolorder(camera_result, c(
  "analysis_id", "signature_id", "inference_role", "formal_axis",
  "formal_n_genes", "formal_direction", "formal_p", "formal_bh_fdr",
  "camera_inter_gene_correlation", "multiplicity_family_n",
  "core_contributing_genes", "context_status", "context_n_genes",
  "context_direction", "context_p_descriptive", "context_bh_fdr_descriptive"
))
fwrite(camera_result, signature_output, sep = "\t", quote = FALSE, na = "NA")
fwrite(rbindlist(core_rows), core_output, sep = "\t", quote = FALSE, na = "NA")

registry <- fromJSON(registry_path, simplifyVector = FALSE)
registry$compact_signature_outputs <- list(
  registry = list(path = signature_registry_output, sha256 = sha256(signature_registry_output)),
  results = list(path = signature_output, sha256 = sha256(signature_output)),
  core_genes = list(path = core_output, sha256 = sha256(core_output)),
  primary_signatures_fdr_lt_0_05 = significant
)
registry$execution_status <- "OVERLAP_AND_SIGNATURES_COMPLETE"
write_json(registry, registry_path, pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(list(status = "PASS", signatures = nrow(camera_result),
                formal_fdr_lt_0_05 = length(significant)), auto_unbox = TRUE), "\n")
