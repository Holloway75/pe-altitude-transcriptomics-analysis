#!/usr/bin/env Rscript

# WARS1 Whole Blood eQTL - PE colocalization (post-selection audit, 2026-08-02).
# Frozen outputs live in
#   results/04_core_candidates/analyses/m4_wars1_post_selection_audit_v1_20260802/colocalization/
# This archived script reads the frozen exposure extract and the frozen M1.1
# GWAS release; reruns write to work/ and must be COMPARED against the frozen
# outputs (see scripts/04_core_candidates/validate_wars1_audit.py), never
# written back into the release.

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

release_dir <- "results/04_core_candidates/analyses/m4_wars1_post_selection_audit_v1_20260802"
eqtl_path <- file.path(release_dir, "colocalization/eqtl_extracts/gtex_v8_whole_blood_WARS1_local.txt")
gwas_path <- "data/processed/01_pe_genetic_map/M1.1_prepare_gwas/FIGSHARE_22680904_v2/metaxcan.tsv.gz"
out_dir <- file.path("work/04_core_candidates/reruns/wars1_coloc_abf", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

probe_bp <- 100800125L
window <- 500000L
region_start <- probe_bp - window
region_end <- probe_bp + window
exp_n <- 670
out_cases <- 16349
out_controls <- 595135
out_n <- out_cases + out_controls

complement <- function(x) {
  chartr("ACGT", "TGCA", x)
}

eqtl <- fread(eqtl_path)
eqtl <- eqtl[
  Chr == 14L & BP >= region_start & BP <= region_end &
    is.finite(b) & is.finite(SE) & SE > 0
]
setorder(eqtl, SNP, p)
eqtl <- eqtl[!duplicated(SNP)]

gwas <- fread(gwas_path)
gwas <- gwas[
  chromosome == 14L & position >= region_start & position <= region_end &
    !is.na(rsid) & rsid != "" &
    is.finite(beta) & is.finite(standard_error) & standard_error > 0 &
    is.finite(effect_allele_frequency)
]
setorder(gwas, rsid, pvalue)
gwas <- gwas[!duplicated(rsid)]

x <- merge(
  eqtl[, .(
    SNP, BP, A1_exp = toupper(A1), A2_exp = toupper(A2),
    beta_exp = b, se_exp = SE, p_exp = p
  )],
  gwas[, .(
    SNP = rsid, BP_out = position,
    EA_out = toupper(effect_allele), OA_out = toupper(non_effect_allele),
    beta_out = beta, se_out = standard_error, p_out = pvalue,
    eaf_out = effect_allele_frequency
  )],
  by = "SNP"
)

x[, match_type := fifelse(
  A1_exp == EA_out & A2_exp == OA_out, "direct",
  fifelse(
    A1_exp == OA_out & A2_exp == EA_out, "swapped",
    fifelse(
      nchar(A1_exp) == 1L & nchar(A2_exp) == 1L &
        A1_exp == complement(EA_out) & A2_exp == complement(OA_out),
      "complement",
      fifelse(
        nchar(A1_exp) == 1L & nchar(A2_exp) == 1L &
          A1_exp == complement(OA_out) & A2_exp == complement(EA_out),
        "complement_swapped", "mismatch"
      )
    )
  )
)]
x <- x[match_type != "mismatch" & BP == BP_out]
x[, beta_out_aligned := fifelse(
  match_type %in% c("swapped", "complement_swapped"), -beta_out, beta_out
)]
x[, maf := pmin(eaf_out, 1 - eaf_out)]
x <- x[is.finite(maf) & maf > 0.01 & maf < 0.5]
setorder(x, BP, SNP)

if (nrow(x) < 50L) {
  stop("Fewer than 50 harmonized variants: ", nrow(x))
}

dataset1 <- list(
  beta = x$beta_exp,
  varbeta = x$se_exp^2,
  snp = x$SNP,
  position = x$BP,
  type = "quant",
  N = exp_n,
  MAF = x$maf
)
dataset2 <- list(
  beta = x$beta_out_aligned,
  varbeta = x$se_out^2,
  snp = x$SNP,
  position = x$BP,
  type = "cc",
  N = out_n,
  s = out_cases / out_n,
  MAF = x$maf
)

prior_grid <- data.table(
  prior_id = c("conservative", "standard", "liberal"),
  p1 = 1e-4,
  p2 = 1e-4,
  p12 = c(1e-6, 1e-5, 1e-4)
)

summaries <- rbindlist(lapply(seq_len(nrow(prior_grid)), function(i) {
  pri <- prior_grid[i]
  fit <- coloc.abf(
    dataset1, dataset2,
    p1 = pri$p1, p2 = pri$p2, p12 = pri$p12
  )
  sm <- fit$summary
  data.table(
    prior_id = pri$prior_id,
    p1 = pri$p1,
    p2 = pri$p2,
    p12 = pri$p12,
    nsnps = as.integer(sm[["nsnps"]]),
    pp_h0 = as.numeric(sm[["PP.H0.abf"]]),
    pp_h1 = as.numeric(sm[["PP.H1.abf"]]),
    pp_h2 = as.numeric(sm[["PP.H2.abf"]]),
    pp_h3 = as.numeric(sm[["PP.H3.abf"]]),
    pp_h4 = as.numeric(sm[["PP.H4.abf"]]),
    h4_given_h3_h4 = as.numeric(sm[["PP.H4.abf"]]) /
      (as.numeric(sm[["PP.H3.abf"]]) + as.numeric(sm[["PP.H4.abf"]]))
  )
}))

regional <- data.table(
  metric = c(
    "region_chr", "region_start", "region_end", "eqtl_rows_in_region",
    "gwas_rows_in_region", "harmonized_maf_gt_0.01", "min_eqtl_p",
    "lead_eqtl_snp", "min_gwas_p", "lead_gwas_snp"
  ),
  value = c(
    "14", region_start, region_end, nrow(eqtl), nrow(gwas), nrow(x),
    format(min(x$p_exp), scientific = TRUE), x$SNP[which.min(x$p_exp)],
    format(min(x$p_out), scientific = TRUE), x$SNP[which.min(x$p_out)]
  )
)

fwrite(x, file.path(out_dir, "wars1_coloc_harmonized.tsv"), sep = "\t")
fwrite(summaries, file.path(out_dir, "wars1_coloc_abf_prior_sensitivity.tsv"), sep = "\t")
fwrite(regional, file.path(out_dir, "wars1_coloc_region_qc.tsv"), sep = "\t")

print(regional)
print(summaries)
