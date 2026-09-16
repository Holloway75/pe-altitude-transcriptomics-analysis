# SNP Matching Rules (Shared Protocol), v1.0

> **Provenance note (2026-09-09).** The original protocol file was lost with the
> removal of the project `docs/` tree (original sha256
> `2b427ada56d586bdf419d9b36a1e51a6235c27d9f3e76258c05549e132f71c13`, as recorded
> in `results/01_pe_genetic_map/analyses/mrjti_native_bonferroni_hg19_v3_20260719/method.json`
> → `matching_protocol`). This document is a reconstruction from (a) the
> `matching_protocol` block of that method.json and (b) the frozen, unchanged
> matching implementations in `scripts/01_pe_genetic_map/` (`prepare_gwas.py`,
> `mrjti_minimal.py`, `build_*_smultixcan_covariance.py`). Declared rules are
> identical; the byte-level document hash differs from the lost original. The
> governing executable artefacts are the frozen scripts and their recorded
> input/output hashes, not this narrative.

## 1. Scope and precedence

This is the single shared SNP-matching specification for GWAS normalization,
GTEx/cis-QTL joins, prediction-model (JTI/S-PrediXcan) alignment, and LD
reference construction in this project. Individual stages must not implement a
divergent matcher. The declared version is `v1.0`.

## 2. Variant identity

- Reference identity is the normalized pair `CHR:POS:REF:ALT` on the target
  build (`GRCh37/hg19` for module 1). Position is 1-based, chromosome without
  `chr` prefix.
- `rsid_policy=required_by_downstream`: canonical, unique `rs[0-9]+` identifiers
  are required in addition to the normalized identity. An rsID never overrides a
  build, position, REF, or ALT conflict.
- Exact normalized `CHR:POS:REF:ALT` keys close every join. Effect/non-effect
  allele ordering is tracked separately from identity and controls the beta
  sign (direct or swapped ordering).

## 3. Reference validation and normalization

- `reference_fasta` = `$HG19_FASTA`
  (`/home/holloway/data/RefGenome/GRCh37/hs37d5/hs37d5.fa`); positive-strand REF
  validation is required before formal matching.
- `coordinate_conversion=none` for registered hg19 inputs. A separately declared
  hg38 intake may use the versioned `rsID_bridge_1KG_hg19` path: canonical rsID,
  same chromosome, and exactly one hg19 1kG `CHR:POS:REF:ALT` candidate, then
  all remaining checks. No implicit liftOver; no rsID-only bridge that bypasses
  allele/FASTA closure.
- `indel_normalization=yes`: left-align and normalize indels against the
  reference.
- `multiallelic_mode=remove_all` (default, incl. the fixed LD background):
  drop every site with more than one ALT. Exception — S-MultiXcan
  common-reference covariance only — protocol mode B
  `multiallelic_mode=exact_alt_match`: require exactly one unique PVAR REF/ALT
  match to the biallelic model (comma-separated ALT and split same-position
  PVAR rows both supported); reject zero or multiple matches.
- Retain normalized biallelic SNPs, indels, and explicitly labelled
  MNV/complex substitutions; exclude symbolic/non-sequence alleles.

## 4. Strand and alleles

- `palindromic_policy=remove`: A/T and C/G single-base SNPs are excluded.
- `strand_policy=allow_nonpalindromic_complement`: non-palindromic variants may
  be complemented to the target strand.
- `effect_flip_policy`: `beta = -beta`, `EAF = 1 - EAF`, `OR = 1/OR` when
  flipping to the target effect allele. Source fields are preserved and the
  alignment status is recorded per variant.
- For indels and MNV/complex substitutions: require normalized full sequence
  strings; no strand-complement rescue.

## 5. QC closure and audit

- Every candidate variant receives one shared-protocol match status.
- Non-unique, palindromic, reference-mismatch, multiallelic, ambiguous, and
  allele-conflict records are excluded with reason-specific counts; every
  exclusion is auditable.
- Final-key QC closure is required: each stage reports matched, flipped, and
  excluded counts by reason, and the joined key set must close against both
  inputs of the join.

## 6. Module-1 stage selections

| Stage | Settings within this protocol |
|---|---|
| GWAS normalization (M1.1) | hg19 source; resolve REF/ALT from source coordinates + row alleles + FASTA; left-align; remove multiallelic; remove palindromic; canonical rsID |
| cis-eQTL (GTEx BESD) | all finite eQTL retained (no significance filter); canonical rsID |
| LD background (MR-JTI) | hg19 EUR TSS ±2 Mb; `indel_normalization=yes`; `multiallelic_mode=remove_all`; MAF ≥ 0.01; missingness ≤ 0.02 |
| S-MultiXcan covariance (M1.5) | mode B `exact_alt_match`; model effect-allele orientation of dosages; ddof=1 ordered covariance cells |
