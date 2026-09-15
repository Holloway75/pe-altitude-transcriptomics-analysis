# Analysis code for "Limited directional convergence between human high-altitude transcriptional responses and genetic susceptibility to pre-eclampsia"

This repository contains the complete data-analysis pipeline that produced the results
of the paper (submitted to *Royal Society Open Science*). The five stage directories
under `scripts/` are verbatim copies of the project scripts, unmodified except that
console-message literals were translated to English on 2026-09-15 (functional reruns
confirmed that all numbers are unchanged; see "Verification" below). A companion
`environment.yml` records the main conda environment. Figure- and table-generation
code for the manuscript itself is not part of this repository.

## Contents

| Directory | Stage |
|---|---|
| `scripts/00_data_preparation/pe_gwas/` | PE GWAS summary-statistic preparation: Figshare VCF ingestion, GCST hg38→hg19 lift-over, rsID indexing, conversion verification |
| `scripts/01_pe_genetic_map/` | PE genetic map: GWAS harmonisation, S-PrediXcan per-tissue prediction, S-MultiXcan covariance and genome-wide runs, MR-JTI selection, downstream assembly and validators |
| `scripts/02_altitude_response/` | Altitude transcriptional response: paired limma differential expression for GSE103927 (primary) and GSE333506 (context), cross-dataset comparison, composition sensitivity, module-3 input freeze, validators |
| `scripts/03_bidirectional_convergence/` | Bidirectional convergence tests: directional pathway tests with empirical nulls, gene-overlap scenarios with donor-level sign-flip tests, continuous-altitude sensitivity, independent numeric verification |
| `scripts/04_core_candidates/` | WARS1 post-selection review: colocalisation (coloc ABF, prior sensitivity, nearby genes), placental single-cell pseudo-bulk audit, QC stress tests, GTEx v11 isoform QTL extraction, audit validator |
| `environment.yml` | Conda specification of the main analysis environment (`deg`) |

60 files (42 Python, 12 R, 5 shell, 1 YAML). Not included: `scripts/04_core_candidates/colocalization/`
and its wrapper `run_colocalization.sh` (a backup colocalisation pipeline with known
defects, retained in the project tree for internal history only — all colocalisation
results in the paper come from the frozen audit scripts in `scripts/04_core_candidates/`
that are included here); the per-module `README.md` working notes (consolidated into
this file); Python bytecode caches; and the manuscript figure/table generation scripts.

## Environments and dependencies

Scripts are executed from the project root with the stage directory as working context
unless a script takes explicit path arguments. Four environments cover the pipeline:

1. **`deg`** (specification: `environment.yml`; tested with R 4.5.3, limma 3.66.0,
   Python 3.12). R side: `limma`, `edgeR`, `data.table`, `Matrix`, `ggplot2`,
   `jsonlite`, `optparse`, `fgsea`. Python side: `pandas`, `numpy`, `scipy`,
   `statsmodels`. Used by: all of `scripts/02_altitude_response` (R),
   `scripts/03_bidirectional_convergence/export_gse103927_limma_prior.R` and
   `run_compact_signatures.R`, and the R scripts in `scripts/04_core_candidates`.
2. **`bio_mr`** (tested with R 4.5.3, coloc 5.2.3). R packages `coloc`, `optparse`,
   `glmnet`, `HDCI` (the latter three are required by the MR-JTI runner
   `mrjti_minimal_runner.R`); external binaries `plink2` (2.0.0-a.6.9LM), `gcta64`
   (1.94.1), `tabix`/`bgzip`. Used by: the colocalisation scripts and the audit
   validator in `scripts/04_core_candidates`, and MR-JTI execution.
3. **`metaxcan`** — the MetaXcan software suite providing `SPrediXcan.py` and
   `SMulTiXcan.py`, invoked by the S-PrediXcan/S-MultiXcan driver scripts in
   `scripts/01_pe_genetic_map` (reference implementation, cited in the manuscript).
4. **Python ≥ 3.9** (tested 3.13.13) — standard library only, except
   `numpy`/`scipy` (null distributions and overlap tests in
   `scripts/03_bidirectional_convergence`) and `pandas`/`pyarrow`
   (`scripts/04_core_candidates/extract_wars1_isoform_qtl.py` reads GTEx parquet).

`scripts/01_pe_genetic_map` additionally resolves tool, reference and data locations
through a project configuration file (`config/paths.env`, loaded by
`scripts/lib/load_config.sh`); every location is overridable through environment
variables, and none of the values is required by the other four stages.

## Input data

Every input is public; sources, accessions and versions are enumerated in the
supplementary methods (electronic supplementary material, methods) and the Data
accessibility statement of the paper. R analysis scripts in
`scripts/02_altitude_response` and `scripts/04_core_candidates` take the local data
root as an argument or through the `PE_DATA_ROOT` environment variable (a directory
containing `GEO/` and `GTEx/` subdirectories with the downloaded public files).
S-PrediXcan/S-MultiXcan/MR-JTI resources (GTEx v8 BESD and covariance files, JTI
tissue models, 1000 Genomes reference panels, dbSNP) are addressed through the
configuration variables described above.

## Output convention

Analysis scripts write rerun outputs only under `work/` and never modify the frozen
result trees under `results/`; a rerun is validated by comparing its outputs against
the frozen release. This convention is enforced by assertions inside the validators.

## Entry points

- Stage 00: `run_convert_pe.sh` / `run_convert_figshare.sh` drive the conversion chain;
  `verify_conversion.sh` (exits non-zero on any mismatch) and `test_conversion.py`
  verify the products.
- Stage 01 execution order: `preflight.py` → `prepare_gwas.py` →
  `run_single_gwas_spredixcan.py` → `build_all_gene_smultixcan_covariance.py` →
  `run_all_gene_smultixcan.py` → `validate_existing_all_gene_smultixcan.py` →
  `assemble_m1_downstream.py` (+ `continue_m1_downstream.py`); MR-JTI via
  `mrjti_minimal.py`. (`run_module01.sh` exposes `status`/`validate`/`publish`
  registration commands; its former `run`/`resume` entry is retired.)
- Stage 02: `analyze_gse103927.R --geo-dir <dir> --hgnc <file> --frozen-dir <dir>
  --output-dir <dir>` (likewise `analyze_gse333506.R`);
  `compare_gse103927_gse333506.R`; `summarize_gse103927_composition_sensitivity.R`;
  `validate_module2_v2_bh.py` for independent recomputation.
- Stage 03: `run_v5.sh` drives the pathway and overlap analyses
  (`run_alt_continuous_threaded.py`, `run_v5_overlap.py`, `validate_v5.py`);
  `verify_v6_numeric_independent.py` performs the independent numeric verification;
  `package_v4_reproduction.py` reproduces the exploratory v4 audit.
- Stage 04: `run_wars1_coloc_abf.R`, `run_nearby_gene_coloc.R`,
  `run_gse173193_wars1_targeted.R`, `qc_stress_gse173193.R`,
  `extract_wars1_isoform_qtl.py`; `validate_wars1_audit.py` checks the frozen audit
  release (7/7 checks).

## Verification

Verified 2026-09-15 on the project tree from which this repository was copied:

- **Syntax**: all 60 files parse cleanly (42 Python via `py_compile`, 12 R via
  `parse()`, 5 shell via `bash -n`).
- **Functional reruns (exact match to the frozen published results)**:
  `analyze_gse103927.R` reproduced the primary analysis — 21 paired donors, 18 224
  genes, 3 183 differentially expressed genes at FDR < 0.05 (1 405 up, 1 778 down);
  WARS1 log~2~FC 0.6421, *P* = 4.36 × 10^-8^, FDR = 1.15 × 10^-5^.
  `analyze_gse333506.R` reproduced the context analysis — 7 paired donors, 12 101
  genes, 2 103 differentially expressed genes.
- **Independent validators**: `validate_module2_v2_bh.py` PASS (7/7 checks,
  Benjamini–Hochberg recomputation); `verify_v6_numeric_independent.py` PASS (8/8
  checks); `validate_wars1_audit.py` PASS (7/7 checks); `preflight.py` PASS;
  `test_mrjti_vendor.py` PASS (4/4 equivalence checks against the upstream MR-JTI
  implementation).
- **Syntax-only**: `compare_gse103927_gse333506.R` (its MSigDB gene-set input file is
  no longer held locally in the exact version used) and `package_v4_reproduction.py`
  (a compute-heavy Monte-Carlo audit whose conclusion-level reproduction is recorded
  in the frozen release) were verified by parsing only.

## Third-party software

The vendored third-party tools used by stage 01 (MR-JTI, MetaXcan, SMR) are not
redistributed in this repository; they are cited in the manuscript references and are
publicly available from their original repositories. No other third-party software is
required beyond the environments listed above.
