#!/usr/bin/env Rscript

# Nearby-gene colocalization comparison (post-selection audit, 2026-08-02).
# Frozen outputs live in
#   results/04_core_candidates/analyses/m4_wars1_post_selection_audit_v1_20260802/colocalization/
# This archived script reads the frozen exposure extracts and the frozen M1.1
# GWAS release; reruns write to work/ and must be COMPARED against the frozen
# outputs, never written back into the release.

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

release_dir <- "results/04_core_candidates/analyses/m4_wars1_post_selection_audit_v1_20260802"
extract_dir <- file.path(release_dir, "colocalization/eqtl_extracts")
out_dir <- file.path("work/04_core_candidates/reruns/nearby_gene_coloc", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gwas <- fread("data/processed/01_pe_genetic_map/M1.1_prepare_gwas/FIGSHARE_22680904_v2/metaxcan.tsv.gz")
gwas <- gwas[chromosome == 14L & position >= 100300125L & position <= 101300125L &
  !is.na(rsid) & is.finite(beta) & is.finite(standard_error) & standard_error > 0 &
  is.finite(effect_allele_frequency)]
setorder(gwas, rsid, pvalue)
gwas <- gwas[!duplicated(rsid)]

comp <- function(x) chartr("ACGT", "TGCA", x)
probes <- data.table(
  probe = c("ENSG00000140105", "ENSG00000197119", "ENSG00000258666", "ENSG00000176473"),
  gene = c("WARS1", "SLC25A29", "RP11-638I2.8", "WDR25"),
  path = file.path(extract_dir, c(
    "gtex_v8_whole_blood_WARS1_local.txt",
    "gtex_v8_whole_blood_ENSG00000197119.txt",
    "gtex_v8_whole_blood_ENSG00000258666.txt",
    "gtex_v8_whole_blood_ENSG00000176473.txt"
  ))
)

res <- rbindlist(lapply(seq_len(nrow(probes)), function(i) {
  e <- fread(probes$path[i])
  e <- e[Chr == 14L & BP >= 100300125L & BP <= 101300125L & is.finite(b) & is.finite(SE) & SE > 0]
  setorder(e, SNP, p)
  e <- e[!duplicated(SNP)]
  x <- merge(
    e[, .(SNP, BP, A1 = toupper(A1), A2 = toupper(A2), b1 = b, se1 = SE, p1 = p)],
    gwas[, .(SNP = rsid, BP2 = position, EA = toupper(effect_allele), OA = toupper(non_effect_allele),
             b2 = beta, se2 = standard_error, p2 = pvalue, eaf = effect_allele_frequency)], by = "SNP"
  )
  x[, mt := fifelse(A1 == EA & A2 == OA, "d",
    fifelse(A1 == OA & A2 == EA, "s",
      fifelse(A1 == comp(EA) & A2 == comp(OA), "c",
        fifelse(A1 == comp(OA) & A2 == comp(EA), "cs", "m"))))]
  x <- x[mt != "m" & BP == BP2]
  x[mt %in% c("s", "cs"), b2 := -b2]
  x[, maf := pmin(eaf, 1 - eaf)]
  x <- x[maf > 0.01 & maf < 0.5]
  d1 <- list(beta=x$b1, varbeta=x$se1^2, snp=x$SNP, position=x$BP,
             type="quant", N=670, MAF=x$maf)
  d2 <- list(beta=x$b2, varbeta=x$se2^2, snp=x$SNP, position=x$BP,
             type="cc", N=611484, s=16349/611484, MAF=x$maf)
  fit <- coloc.abf(d1, d2, p1=1e-4, p2=1e-4, p12=1e-5)
  z <- fit$summary
  data.table(probe=probes$probe[i], gene=probes$gene[i], nsnps=z[["nsnps"]],
    lead_eqtl_snp=x$SNP[which.min(x$p1)], min_eqtl_p=min(x$p1),
    pp_h3=z[["PP.H3.abf"]], pp_h4=z[["PP.H4.abf"]],
    h4_given_h3_h4=z[["PP.H4.abf"]]/(z[["PP.H3.abf"]]+z[["PP.H4.abf"]]))
}))

fwrite(res, file.path(out_dir, "nearby_gene_coloc_standard_prior.tsv"), sep="\t")
print(res)
