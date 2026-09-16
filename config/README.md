# Configuration

Run parameters are centralised in this directory; analysis scripts never embed
personal absolute paths.

- `paths.env` — machine-specific locations: external data root (JTI models,
  MSigDB, MAGMA, reference panels), software (conda environments, binaries) and
  the project-internal `work/`, `data/processed/` and `results/` directories.
  When moving to another machine, either edit this file or export environment
  variables of the same names before running; variables always take precedence
  over the file defaults. `scripts/lib/load_config.sh` loads it (and
  `analysis.env`) with `set -a`, so shell entry points see every value.
- `analysis.env` — shared analysis parameters, GWAS ordering and thresholds;
  also freezes the module-2 formal/context ALT axes and the module-2→3 release
  ID, and marks the archived module-3 v4 analysis `ARCHIVED_NEGATIVE`
  (historical metadata only; it does not license active execution).
- `pe_gwas.tsv` — phenotype, sample size, build and analysis role of the four
  PE GWAS datasets.
- `expression_datasets.tsv` — explicit roles of the formal GSE103927 ALT axis,
  the GSE333506 external-context axis and the archived/retired GEO datasets.
- `module01.json` — module 1 schema, candidate rules, scientific QC,
  cross-tissue status and gene-annotation/LD-block resource paths.
- `module03_v5.json` / `module03_v6.json` — module 3 (bidirectional
  convergence) parameters. The frozen published analysis is
  `m3_exploratory_shared_programs_v6_20260910` (driven by `module03_v6.json`
  through `run_v5_overlap.py` / `run_compact_signatures.R` / `validate_v5.py`);
  `run_v5.sh` documents the earlier v5 entry and reads `module03_v5.json`.

R/Python analysis scripts receive explicit command-line arguments from the
entry points; they do not assume local paths on their own. The main internal
directory variables are `INTERIM_DIR` (`work/`), `PROCESSED_DIR`,
`RESULTS_DIR` and the per-module derived paths.
