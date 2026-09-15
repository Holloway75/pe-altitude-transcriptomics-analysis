#!/usr/bin/env Rscript

# GSE173193 placental single-cell WARS pseudo-bulk audit (post-selection, 2026-08-02).
# Frozen outputs live in
#   results/04_core_candidates/analyses/m4_wars1_post_selection_audit_v1_20260802/single_cell/
# Reruns write to work/ and must be COMPARED against the frozen outputs,
# never written back into the release.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

# Local data root containing GEO/ and GTEx/ subdirectories downloaded from the
# public sources; override with the PE_DATA_ROOT environment variable.
data_dir <- file.path(Sys.getenv("PE_DATA_ROOT", unset = "data"),
                      "GEO", "GSE173193_minimal_PE_control")
out_dir <- file.path("work/04_core_candidates/reruns/gse173193_wars1_targeted", format(Sys.Date(), "%Y%m%d"))
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

read_sample <- function(accession, sample, condition) {
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
  keep <- detected >= 200 & lib >= 500 & is.finite(mt_fraction) & mt_fraction < 0.20
  mat <- mat[, keep, drop = FALSE]
  lib <- lib[keep]
  detected <- detected[keep]
  mt_fraction <- mt_fraction[keep]

  all_markers <- unique(c(unlist(markers), unlist(troph_markers), "WARS"))
  present <- intersect(all_markers, rownames(mat))
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

  wars <- as.numeric(mat["WARS", ])
  data.table(
    accession, sample, condition, barcode = colnames(mat), cell_type = broad,
    library_size = lib, detected_genes = detected, mt_fraction,
    wars_count = wars, wars_log1p_cpm = log1p(wars / lib * 1e4)
  )
}

cells <- rbindlist(lapply(seq_len(nrow(samples)), function(i) {
  read_sample(samples$accession[i], samples$sample[i], samples$condition[i])
}))

sample_celltype <- cells[, .(
  n_cells = .N,
  wars_detected_fraction = mean(wars_count > 0),
  wars_mean_log1p_cpm = mean(wars_log1p_cpm),
  wars_pseudobulk_cpm = sum(wars_count) / sum(library_size) * 1e6
), by = .(condition, sample, accession, cell_type)]

overall_sample <- cells[, .(
  n_cells = .N,
  wars_detected_fraction = mean(wars_count > 0),
  wars_mean_log1p_cpm = mean(wars_log1p_cpm),
  wars_pseudobulk_cpm = sum(wars_count) / sum(library_size) * 1e6
), by = .(condition, sample, accession)]
overall_direction <- overall_sample[, .(
  n_donors = .N,
  mean_pseudobulk_cpm = mean(wars_pseudobulk_cpm),
  donor_values = paste(sample, sprintf("%.3f", wars_pseudobulk_cpm), sep = ":", collapse = ";")
), by = condition]
overall_direction <- dcast(overall_direction, . ~ condition,
  value.var = c("n_donors", "mean_pseudobulk_cpm", "donor_values"))
overall_direction[, log2_pe_vs_control := log2((mean_pseudobulk_cpm_PE + 0.5) /
                                               (mean_pseudobulk_cpm_Control + 0.5))]

direction <- sample_celltype[n_cells >= 20, .(
  n_donors = uniqueN(sample),
  min_cells_per_donor = min(n_cells),
  mean_pseudobulk_cpm = mean(wars_pseudobulk_cpm),
  donor_values = paste(sample, sprintf("%.3f", wars_pseudobulk_cpm), sep = ":", collapse = ";")
), by = .(condition, cell_type)]
direction <- dcast(direction, cell_type ~ condition,
  value.var = c("n_donors", "min_cells_per_donor", "mean_pseudobulk_cpm", "donor_values"))
direction[, log2_pe_vs_control := log2((mean_pseudobulk_cpm_PE + 0.5) /
                                       (mean_pseudobulk_cpm_Control + 0.5))]
direction[, direction := fifelse(log2_pe_vs_control > 0, "PE_higher",
                                 fifelse(log2_pe_vs_control < 0, "PE_lower", "equal"))]

qc <- cells[, .(
  cells_after_qc = .N,
  median_library_size = median(library_size),
  median_detected_genes = median(detected_genes),
  median_mt_fraction = median(mt_fraction)
), by = .(condition, sample, accession)]

fwrite(qc, file.path(out_dir, "gse173193_sample_qc.tsv"), sep = "\t")
fwrite(overall_sample, file.path(out_dir, "gse173193_wars1_overall_by_sample.tsv"), sep = "\t")
fwrite(overall_direction, file.path(out_dir, "gse173193_wars1_overall_direction.tsv"), sep = "\t")
fwrite(sample_celltype, file.path(out_dir, "gse173193_wars1_by_sample_celltype.tsv"), sep = "\t")
fwrite(direction, file.path(out_dir, "gse173193_wars1_direction.tsv"), sep = "\t")
fwrite(cells[, .N, by = .(condition, sample, cell_type)],
       file.path(out_dir, "gse173193_marker_celltype_counts.tsv"), sep = "\t")

print(qc)
print(overall_sample)
print(overall_direction)
print(direction)
