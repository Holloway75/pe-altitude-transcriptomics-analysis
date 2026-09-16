# Analysis code for "Limited directional convergence between human high-altitude transcriptional responses and genetic susceptibility to pre-eclampsia"

This repository contains the complete data-analysis pipeline that produced the results
of the paper (submitted to *Royal Society Open Science*). The stage directories under
`scripts/`, the shared `config/` tree and `scripts/lib/`, and the protocol note under
`docs/` are verbatim copies of the project files. Relative to the original 2026-09-15
snapshot (console-message literals translated to English; see "Verification" below),
a 2026-09-16 repair release closed the reproducibility gaps found by a full-pipeline
rerun audit — see "Repository fixes (2026-09-16)" for the exact changes. A companion
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
| `config/` | Shared run configuration: machine paths (`paths.env`), analysis parameters (`analysis.env`), dataset registries (`pe_gwas.tsv`, `expression_datasets.tsv`), module parameters (`module01.json`, `module03_v5.json`, `module03_v6.json`) |
| `scripts/lib/` | `load_config.sh`, the shell configuration loader sourced by the stage entry points |
| `docs/protocol/` | `SNP_MATCHING_RULES.md`, the shared SNP-matching specification hashed into the MR-JTI method records |
| `environment.yml` | Conda specification of the main analysis environment (`deg`), pinned to the tested versions |

73 files (43 Python, 12 R, 6 shell, 7 configuration, 3 Markdown, 2 YAML). A second,
unpinned `environment.yml` inside `scripts/02_altitude_response/` records the fuller
historical stage-02 workspace (including GEO download-era packages) and is retained
from the original snapshot. Not included: `scripts/04_core_candidates/colocalization/`
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

### Configuration (`config/`, `scripts/lib/load_config.sh`)

Machine-specific locations resolve through `config/paths.env`, which is sourced
(via `scripts/lib/load_config.sh`) by the stage-00 shell drivers, the stage-03
entry script, and by the stage-01 Python drivers through a `bash -c "set -a;
source scripts/lib/load_config.sh"` environment snapshot. `config/analysis.env`
fixes shared analysis parameters and GWAS ordering; `config/pe_gwas.tsv` and
`config/expression_datasets.tsv` are the dataset registries; `config/module01.json`
and `config/module03_v5.json`/`module03_v6.json` carry the module parameters
hashed into the method records. Every location default in `paths.env` can be
overridden by an environment variable of the same name (for example
`EXTERNAL_DATA_ROOT=/mnt/data`), so the files publish the reference layout
without binding a user to it. See `config/README.md`.

## Input data

Every input is public; sources, accessions and versions are enumerated in the
supplementary methods (electronic supplementary material, methods) and the Data
accessibility statement of the paper. Concretely:

- **Expression data (stages 02–04)**: R analysis scripts in
  `scripts/02_altitude_response` and `scripts/04_core_candidates` take the local
  data root as an argument or through the `PE_DATA_ROOT` environment variable (a
  directory containing `GEO/` and `GTEx/` subdirectories with the downloaded
  public files).
- **GWAS and genetic-map resources (stages 00–01)**: addressed through
  `config/paths.env`. The reference external-data layout is
  `GWAS_summary/05_outcomes/PE/` (the four PE GWAS, Figshare + three GCST),
  `RefGenome/GRCh37/hs37d5/hs37d5.fa` and `RefGenome/GRCh38/` (reference FASTAs
  with `.fai` indices), `dbSNP/hg38_build157/`, `JTI_model/` (49 paired
  `JTI_<tissue>.db` + `.txt.gz` covariance files), `1kGenome/` (hg19/hg38 PLINK
  binary EUR panels and hg19 VCFs) and `MSigDB/2024.1.Hs/`.
- **Third-party tools (stage 01)** are not redistributed; place them under
  `tools/` in the project root so that the hardcoded relative paths resolve:
  `tools/MetaXcan/software/SPrediXcan.py` and `SMulTiXcan.py`, `tools/MR-JTI/mr/MR-JTI.r`
  (equivalent to upstream gamzonlab/MR-JTI master; see
  `scripts/01_pe_genetic_map/test_mrjti_vendor.py`), and `tools/smr/smr`.
  `scripts/01_pe_genetic_map/preflight.py` checks for exactly these files.
- **Gene nomenclature snapshot (stage 02)**: the analyses were run against
  `hgnc_complete_set_current_20260809.tsv`
  (sha256 `faaeb6ae1e2a596be658b5f23ee44937c9c5379fa37d1f5ace54b94b16962b1c`).
  HGNC serves only the current snapshot, so the exact file cannot be re-downloaded
  from the upstream site; the recorded hash lets an archived copy be verified.
- **Gene sets (stages 02–03)**: MSigDB 2024.1.Hs symbol GMTs —
  `h.all.v2024.1.Hs.symbols.gmt`
  (sha256 `ee2463540042078bfa3f67828e1e223bb354446d9fbb4d22845866835ba5c772`)
  and `c2.cp.reactome.v2024.1.Hs.symbols.gmt`
  (sha256 `9ea1b5e656597daf423e41c5ebcaa9892bfedf3292fff768605d7b0d5e5e9703`) —
  both still downloadable from the GSEA/MSigDB archive of past releases.

## Output convention

Analysis scripts write rerun outputs only under `work/` and never modify the frozen
result trees under `results/`; a rerun is validated by comparing its outputs against
the frozen release. This convention is enforced by assertions inside the validators.

## Entry points

- Stage 00: `run_convert_pe.sh` / `run_convert_figshare.sh` drive the conversion chain
  (both source `scripts/lib/load_config.sh`, so `config/` must be present);
  `verify_conversion.sh` (exits non-zero on any mismatch) and `test_conversion.py`
  verify the products.
- Stage 01 execution order: `preflight.py` → `prepare_gwas.py` →
  `run_single_gwas_spredixcan.py` → `validate_spredixcan_release.py` (writes the
  per-dataset M1.3 `validation-report.json`) → `build_all_gene_smultixcan_covariance.py`
  → `run_all_gene_smultixcan.py` → `validate_existing_all_gene_smultixcan.py` →
  `assemble_m1_downstream.py` (+ `continue_m1_downstream.py`). The published MR-JTI
  chain is `run_global_bh_mrjti.py` (one `--source DATASET=table` per GWAS) →
  `validate_global_bh_mrjti.py <run-dir> --recompute-source-bh` →
  `reinfer_mrjti_native_bonferroni.py --source <run-dir> --output <v3-dir>` →
  `validate_global_bh_mrjti.py <v3-dir>`. (`run_module01.sh` exposes `status`/
  `validate`/`publish` registration commands; its former `run`/`resume` entry is
  retired, as is the integrated `execute.py` run; `mrjti_minimal.py` is the
  single-task MR-JTI entry.)
- Stage 02: `analyze_gse103927.R --geo-dir <dir> --hgnc <file> --frozen-dir <dir>
  --output-dir <dir>` (likewise `analyze_gse333506.R`);
  `compare_gse103927_gse333506.R`; `summarize_gse103927_composition_sensitivity.R`;
  `validate_module2_v2_bh.py` for independent recomputation. (`validate_module2.py`/
  `seal_module2.py` are legacy retired-schema tools; they read
  `config/expression_datasets.tsv`.)
- Stage 03: the frozen published analysis is v6
  (`m3_exploratory_shared_programs_v6_20260910`), driven by
  `run_v5_overlap.py --config config/module03_v6.json` →
  `run_compact_signatures.R --config config/module03_v6.json` →
  `validate_v5.py`; `export_gse103927_limma_prior.R` regenerates the module-2
  prior input; `verify_v6_numeric_independent.py` performs the independent numeric
  verification; `package_v4_reproduction.py` reproduces the exploratory v4 audit.
  `run_v5.sh` documents the earlier v5 entry (it additionally requires the archived
  `run_alt_continuous.py`, which is not part of this repository). A full stage-03
  rerun additionally needs the module-2 release directory and the v4 gene-universe/
  mapping registry from the corresponding frozen project releases.
- Stage 04: `run_wars1_coloc_abf.R`, `run_nearby_gene_coloc.R`,
  `run_gse173193_wars1_targeted.R`, `qc_stress_gse173193.R`,
  `extract_wars1_isoform_qtl.py`; `validate_wars1_audit.py` checks the frozen audit
  release (7/7 checks; the coloc recomputation check must be run with the `bio_mr`
  Rscript, which provides the `coloc` package).

## Repository fixes (2026-09-16)

A from-scratch rerun of the whole pipeline from this repository (isolated sandbox,
all five stages, every output compared against the frozen release — see
"Verification") found no scientific discrepancies, but three defects that stopped a
third-party rerun before it could reach the published numbers, plus several
documentation gaps. This release fixes them:

1. **`config/`, `scripts/lib/` and `docs/protocol/SNP_MATCHING_RULES.md` are now
   included.** They are hard requirements of the stage-00/01/03 entry points (the
   drivers source `scripts/lib/load_config.sh`, which fails without
   `paths.env`/`analysis.env`; the MR-JTI runner and reinferencer hash
   `docs/protocol/SNP_MATCHING_RULES.md` into their method records), but were not
   part of the original snapshot. The published files are the project files with
   only two comment lines in `paths.env` translated from Chinese.
2. **`scripts/01_pe_genetic_map/validate_spredixcan_release.py` is new** (it also
   exists in the project tree). The frozen M1.3 releases each contain a
   `validation-report.json` that `reinfer_mrjti_native_bonferroni.py` hard-requires,
   but no script in the original snapshot could produce it. The new script
   recomputes every reported quantity from the M1.3 artifacts (including a full
   Benjamini–Hochberg recomputation over the association table) and writes the
   report in the frozen schema; run it after `run_single_gwas_spredixcan.py`.
   Verified against the four frozen reports: field-for-field identical except the
   `method.json` hash, which necessarily changes with each run's paths.
3. **`reinfer_mrjti_native_bonferroni.py` accepted only the legacy status string**
   (`success_inference_descriptive`), so on the output of the current
   `run_global_bh_mrjti.py` (which inline-enriches rows as
   `success_inference_bonferroni_selection_dependent`) it nulled
   `mrjti_supported`/`ci_significance` and the subsequent validator failed. It now
   recomputes the native-Bonferroni CIs from the frozen bootstrap draws for either
   status. Verified: reinfer on a fresh runner output reproduces the frozen v3
   table with zero differing cells (2 335 rows, 78 supported) and the final
   `validate_global_bh_mrjti.py` passes (24 053 artifacts, 0 errors).
4. **`scripts/02_altitude_response/test_validate_module2.py`** skipped one test
   that reads an internal project reference file when that file is absent, so the
   suite is self-contained (35 tests: 34 pass, 1 documented skip in this copy).
5. **`environment.yml` is now actually present** (the original README referenced
   it, but the file was missing from the snapshot).
6. **README**: corrected the configuration claim (all shell entry points need
   `config/`, not only stage 01), documented the data-root layout, the `tools/`
   layout, the exact HGNC/MSigDB snapshots with hashes, the MR-JTI chain and the
   v6 stage-03 entry, and added the known-differences list below.

## Known differences from the frozen release

A rerun with this code reproduces every scientific value, but not every byte. The
known, explained differences are:

- `gene_execution_audit.tsv` (M1.5): 1 461 rows with no model output carry the
  label `outcome_unavailable` where the frozen release wrote
  `exception_or_missing_output`; the success/degenerate/underflow counts are
  identical and no downstream number changes.
- MR-JTI CI text precision: the runner writes 15-digit R text for
  `ci_low`/`ci_high`; the frozen v3 layer and the reinferencer use `.17g`. The
  values are the same doubles; running the (fixed) reinferencer restores the
  frozen formatting exactly.
- Container timestamps/encodings: python-gzip headers carry an mtime (the frozen
  M1.1 `metaxcan.tsv.gz` are content-identical, not byte-identical), `.npz` zip
  container bytes differ by timestamps with equal arrays, and current htslib
  writes `.tbi` indices 20–30 KB larger than the archived ones with identical
  VCF content.
- Assembly provenance columns (`pe_cross_tissue.tsv`, `pe_candidate_locus_map.tsv`):
  input-hash columns embed run paths/timestamps and differ; all value columns are
  equal.
- Bookkeeping: the M1.1 `method.json` `datasets.*.input_sha256` values recorded by
  the frozen run do not match any hash derivable from the shipped inputs (the
  output-side hash chain is complete); the archived module-4 validation report
  records 6/7 because its coloc recomputation was run with an Rscript lacking the
  `coloc` package at archive time — with the specified `bio_mr` Rscript all 7/7
  checks pass (the frozen directory itself is unchanged).

## Verification

Verified 2026-09-15 on the project tree from which this repository was copied:

- **Syntax**: all files parse cleanly (Python via `py_compile`, R via `parse()`,
  shell via `bash -n`).
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

Verified 2026-09-16 by a full-pipeline rerun from this repository in an isolated
sandbox (fresh clone plus documented bridges for the non-redistributed inputs;
rerun outputs compared against the frozen releases):

- All five stages reproduced every scientific value: stage-00 conversions
  (three GCST hg19 products byte-identical; Figshare product content-identical,
  verification PASS), stage-01 M1.1/M1.3/M1.5 tables (byte-identical core tables),
  the M1.7–M1.9 assembly (row counts equal; `pe_whole_blood.tsv` byte-identical
  with the registry-recorded sha256), the MR-JTI batch (2 335 tasks; candidate
  manifest byte-identical; all statistics columns equal), stage-02 limma analyses
  (DE tables byte-identical), the stage-03 v6 chain (five result tables
  byte-identical; 100 000-draw sign-flip matrix arrays equal; validators 15/15 and
  8/8), composition sensitivity (6/6 files byte-identical) and the stage-04 WARS1
  audit (colocalisation, targeted, QC-stress and isoform tables byte-identical;
  validator 7/7).
- The fix verifications listed under "Repository fixes (2026-09-16)".
- One environment note: the Figshare conversion performs per-variant FASTA seeks;
  on cold page cache over slow shared storage this can be extremely slow. Copying
  the reference FASTA to local fast storage before the run restores the expected
  throughput (about 5 minutes for the 6.6 GB VCF).

## Third-party software

The vendored third-party tools used by stage 01 (MR-JTI, MetaXcan, SMR) are not
redistributed in this repository; they are cited in the manuscript references and are
publicly available from their original repositories. No other third-party software is
required beyond the environments listed above.
