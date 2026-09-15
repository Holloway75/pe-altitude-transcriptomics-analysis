#!/usr/bin/env Rscript

# QC stress tests for the GSE173193 WARS pseudo-bulk audit (Section 3.5).
# Read-only with respect to the frozen release; all outputs written under work/.
# Questions addressed:
#  A. QC-threshold sensitivity of the overall PE-vs-control direction
#  B. Composition confounding: direct-standardized overall WARS CPM
#  C. Library-depth confounding at the donor level
#  D. Exact donor-level permutation test (all 6 label allocations)
#  E. Cell-type direction robustness across QC thresholds
#  F. Alternative estimand: mean log1p CPM per donor

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

# Local data root containing GEO/ and GTEx/ subdirectories downloaded from the
# public sources; override with the PE_DATA_ROOT environment variable.
data_dir <- file.path(Sys.getenv("PE_DATA_ROOT", unset = "data"),
                      "GEO", "GSE173193_minimal_PE_control")
out_dir <- "work/04_core_candidates/qc_stress"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- data.table(
  accession = c("GSM5261695", "GSM5261696", "GSM5261699", "GSM5261700"),
  sample = c("C1", "C2", "P1", "P2"),
  condition = c("Control", "Control", "PE", "PE")
)

markers <- list(
  Trophoblast = c("KRT7", "KRT8", "KRT18", "GATA3", "TFAP2A", "TFAP2C", "EPCAM"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG"),
  Myeloid = c("PTPRC", "LST1", "TYROBP", "FCER1G", "CD68", "C1QC"),
  Lymphoid = c("PTPRC", "CD3D", "CD3E", "NKG7", "GNLY"),
  Stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  Erythroid = c("HBB", "HBA1", "HBA2", "GYPA")
)
troph_markers <- list(
  CTB = c("PAGE4", "TP63", "EGFR", "PEG10"),
  EVT = c("HLA-G", "MMP2", "ITGA5", "FN1"),
  SCT = c("CGB3", "CSH1", "PSG1", "SDC1", "ERVW-1")
)

read_sample_raw <- function(accession, sample, condition) {
  prefix <- file.path(data_dir, paste0(accession, "_", sample))
  feature_path <- paste0(prefix, "_features.tsv.gz")
  if (!file.exists(feature_path)) {
    feature_path <- file.path(data_dir, "GSM5261695_C1_features.tsv.gz")
  }
  features <- fread(feature_path, header = FALSE)
  mat <- readMM(gzfile(paste0(prefix, "_matrix.mtx.gz")))
  barcode_path <- paste0(prefix, "_barcodes.tsv.gz")
  barcodes <- if (file.exists(barcode_path)) {
    fread(barcode_path, header = FALSE)$V1
  } else {
    sprintf("cell_%06d", seq_len(ncol(mat)))
  }
  genes <- make.unique(features$V2)
  rownames(mat) <- genes
  colnames(mat) <- barcodes
  lib <- Matrix::colSums(mat)
  detected <- Matrix::colSums(mat > 0)
  mt_idx <- grep("^MT-", genes)
  mt_fraction <- if (length(mt_idx)) Matrix::colSums(mat[mt_idx, , drop = FALSE]) / lib else rep(0, ncol(mat))
  all_markers <- unique(c(unlist(markers), unlist(troph_markers), "WARS"))
  list(mat = mat, genes = genes, lib = lib, detected = detected,
       mt_fraction = mt_fraction, sample = sample, condition = condition,
       accession = accession, all_markers = all_markers)
}

raws <- lapply(seq_len(nrow(samples)), function(i) {
  read_sample_raw(samples$accession[i], samples$sample[i], samples$condition[i])
})

annotate <- function(raw, keep) {
  mat <- raw$mat[, keep, drop = FALSE]
  lib <- raw$lib[keep]
  present <- intersect(raw$all_markers, rownames(mat))
  norm <- log1p(t(t(mat[present, , drop = FALSE]) / lib) * 1e4)
  score <- sapply(markers, function(gs) {
    gs <- intersect(gs, rownames(norm))
    Matrix::colMeans(norm[gs, , drop = FALSE])
  })
  broad_idx <- max.col(score, ties.method = "first")
  broad <- colnames(score)[broad_idx]
  broad[apply(score, 1, max) <= 0.10] <- "Unassigned"
  trop_cells <- which(broad == "Trophoblast")
  if (length(trop_cells)) {
    ts <- sapply(troph_markers, function(gs) {
      gs <- intersect(gs, rownames(norm))
      Matrix::colMeans(norm[gs, trop_cells, drop = FALSE])
    })
    subtype_idx <- max.col(ts, ties.method = "first")
    subtype <- colnames(ts)[subtype_idx]
    subtype[apply(ts, 1, max) <= 0.10] <- "Other"
    broad[trop_cells] <- paste0("Trophoblast_", subtype)
  }
  broad
}

run_threshold <- function(min_genes, min_lib, max_mt) {
  rbindlist(lapply(raws, function(raw) {
    keep <- raw$detected >= min_genes & raw$lib >= min_lib &
      is.finite(raw$mt_fraction) & raw$mt_fraction < max_mt
    ct <- annotate(raw, keep)
    wars <- as.numeric(raw$mat["WARS", keep])
    lib <- raw$lib[keep]
    data.table(
      sample = raw$sample, condition = raw$condition,
      cell_type = ct, wars_count = wars, library_size = lib
    )
  }))
}

log2_pe_ctrl <- function(dt) {
  smp <- dt[, .(cpm = sum(wars_count) / sum(library_size) * 1e6),
            by = .(condition, sample)]
  m <- smp[, mean(cpm), by = condition]
  v <- setNames(m$V1, m$condition)
  if (is.na(v["PE"]) || is.na(v["Control"])) return(NA_real_)
  log2((v["PE"] + 0.5) / (v["Control"] + 0.5))
}

# ---- A. QC threshold grid ------------------------------------------------
# 2026-09-10 (audit A-3 / Minor-6): the grid originally tightened only
# (200/500/0.20 -> stricter). Relaxed levels (100 genes, 250 UMIs, 30% MT)
# are added so both directions of the stress test are covered; the frozen
# release thresholds remain 200 genes / 500 UMIs / 20% MT (baseline row).
grid <- CJ(
  min_genes = c(100, 200, 500, 1000),
  min_lib = c(250, 500, 1000, 2500),
  max_mt = c(0.10, 0.20, 0.30)
)
grid[, log2_pe_vs_control := NA_real_]
grid[, n_cells := NA_integer_]
for (i in seq_len(nrow(grid))) {
  dt <- run_threshold(grid$min_genes[i], grid$min_lib[i], grid$max_mt[i])
  grid[i, log2_pe_vs_control := log2_pe_ctrl(dt)]
  grid[i, n_cells := nrow(dt)]
}
fwrite(grid, file.path(out_dir, "qc_threshold_grid.tsv"), sep = "\t")
print(grid)
# Baseline row must reproduce the frozen section 3.5 overall direction
# (log2 PE/control = -1.17 from run_gse173193_wars1_targeted.R).
baseline_row <- grid[min_genes == 200 & min_lib == 500 & max_mt == 0.20]
stopifnot(nrow(baseline_row) == 1)
stopifnot(abs(baseline_row$log2_pe_vs_control - (-1.17)) < 0.005)

# ---- B. Composition standardization (frozen thresholds) ------------------
base <- run_threshold(200, 500, 0.20)
smp_ct <- base[, .(cpm = sum(wars_count) / sum(library_size) * 1e6,
                   n = .N), by = .(sample, condition, cell_type)]
# reference composition = mean across donors of each cell type's share
share <- base[, .(.N), by = .(sample, cell_type)]
share[, frac := N / sum(N), by = sample]
ref_comp <- share[, .(ref_frac = mean(frac)), by = cell_type]
# standardized CPM per donor = sum over ct of ref_frac * donor_ct_cpm (where present)
std <- merge(smp_ct, ref_comp, by = "cell_type")
std_donor <- std[, .(std_cpm = sum(ref_frac * cpm)), by = .(sample, condition)]
std_group <- std_donor[, .(mean_std_cpm = mean(std_cpm)), by = condition]
std_res <- dcast(std_group, . ~ condition, value.var = "mean_std_cpm")
std_res[, log2_pe_vs_control := log2((PE + 0.5) / (Control + 0.5))]
fwrite(std_donor, file.path(out_dir, "composition_std_donor.tsv"), sep = "\t")
fwrite(std_res, file.path(out_dir, "composition_std_overall.tsv"), sep = "\t")
print(std_donor); print(std_res)

# ---- C. donor-level depth diagnostics ------------------------------------
donor <- base[, .(n_cells = .N,
                  total_lib = sum(library_size),
                  median_lib = median(library_size),
                  wars_cpm = sum(wars_count) / sum(library_size) * 1e6,
                  wars_detected_frac = mean(wars_count > 0)),
              by = .(sample, condition)]
print(donor)
fwrite(donor, file.path(out_dir, "donor_depth_diagnostics.tsv"), sep = "\t")

# ---- D. exact donor permutation (all 6 allocations) ----------------------
vals <- setNames(donor$wars_cpm, donor$sample)
cond <- setNames(donor$condition, donor$sample)
allocs <- combn(names(vals), 2, simplify = FALSE)
perm <- rbindlist(lapply(allocs, function(pe_set) {
  pe_v <- vals[pe_set]; ctrl_v <- vals[setdiff(names(vals), pe_set)]
  data.table(
    pe_donors = paste(pe_set, collapse = "+"),
    log2_ratio = log2((mean(pe_v) + 0.5) / (mean(ctrl_v) + 0.5)),
    matches_observed = setequal(pe_set, names(cond[cond == "PE"]))
  )
}))
setorder(perm, log2_ratio)
fwrite(perm, file.path(out_dir, "exact_donor_permutation.tsv"), sep = "\t")
print(perm)
obs <- perm[matches_observed == TRUE, log2_ratio]
cat(sprintf("two-sided exact p = %.4f (6 allocations)\n",
            mean(abs(perm$log2_ratio) >= abs(obs))))

# ---- E. cell-type direction across thresholds ----------------------------
ct_res <- list()
for (i in seq_len(nrow(grid))) {
  dt <- run_threshold(grid$min_genes[i], grid$min_lib[i], grid$max_mt[i])
  sct <- dt[, .(n_cells = .N, cpm = sum(wars_count) / sum(library_size) * 1e6),
            by = .(condition, sample, cell_type)]
  d <- sct[n_cells >= 20, .(n_donors = uniqueN(sample), cpm = mean(cpm)),
           by = .(condition, cell_type)]
  dc <- dcast(d, cell_type ~ condition, value.var = c("n_donors", "cpm"))
  dc[, log2 := log2((cpm_PE + 0.5) / (cpm_Control + 0.5))]
  dc[, grid := sprintf("g%d_%d_%.2f", grid$min_genes[i], grid$min_lib[i], grid$max_mt[i])]
  ct_res[[i]] <- dc
}
ct_all <- rbindlist(ct_res)
fwrite(ct_all, file.path(out_dir, "celltype_direction_threshold_grid.tsv"), sep = "\t")
key_cts <- c("Trophoblast_CTB", "Trophoblast_EVT", "Endothelial", "Myeloid")
print(ct_all[cell_type %in% key_cts,
             .(grid, cell_type, n_donors_Control, n_donors_PE, log2)])

# ---- F. alternative estimand: mean log1p CPM per donor -------------------
alt <- rbindlist(lapply(raws, function(raw) {
  keep <- raw$detected >= 200 & raw$lib >= 500 &
    is.finite(raw$mt_fraction) & raw$mt_fraction < 0.20
  wars <- as.numeric(raw$mat["WARS", keep])
  lib <- raw$lib[keep]
  data.table(sample = raw$sample, condition = raw$condition,
             mean_log1p_cpm = mean(log1p(wars / lib * 1e4)),
             detected_frac = mean(wars > 0))
}))
altg <- alt[, lapply(.SD, mean), by = condition,
            .SDcols = c("mean_log1p_cpm", "detected_frac")]
print(altg)
fwrite(alt, file.path(out_dir, "alternative_estimand_donor.tsv"), sep = "\t")
fwrite(altg, file.path(out_dir, "alternative_estimand_group.tsv"), sep = "\t")

cat("\nSTRESS TESTS COMPLETE\n")
