#!/usr/bin/env python3
"""Build a single-gene S-MultiXcan SNP-dosage covariance from hg19 PGEN.

The reference PGEN may contain multiallelic records.  Model variants are
matched to one exact REF/ALT allele in the PVAR, then allele-specific dosages
are read from a small PLINK2-exported VCF.  Covariance is calculated on model
effect-allele dosage with ddof=1, matching MetaXcan's CovarianceBuilder.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import re
import sqlite3
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np

RSID = re.compile(r"^rs[0-9]+$")
DNA = re.compile(r"^[ACGT]+$")
COMP = str.maketrans("ACGT", "TGCA")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def complement(base: str) -> str:
    return base.translate(COMP)


def is_palindromic(a1: str, a2: str) -> bool:
    return len(a1) == len(a2) == 1 and complement(a1) == a2


def match_model_alleles(model_ref: str, model_eff: str, pvar_ref: str, pvar_alt: str):
    """Return shared-protocol status and PVAR allele counted by the model."""
    model_ref, model_eff = model_ref.upper(), model_eff.upper()
    pvar_ref, pvar_alt = pvar_ref.upper(), pvar_alt.upper()
    if not all(DNA.fullmatch(x) for x in (model_ref, model_eff, pvar_ref, pvar_alt)):
        return None
    if model_ref == model_eff or pvar_ref == pvar_alt:
        return None
    if is_palindromic(model_ref, model_eff):
        return ("palindrome_removed", None)
    if (model_ref, model_eff) == (pvar_ref, pvar_alt):
        return ("exact", pvar_alt)
    if (model_ref, model_eff) == (pvar_alt, pvar_ref):
        return ("flipped", pvar_ref)
    if len(model_ref) == len(model_eff) == len(pvar_ref) == len(pvar_alt) == 1:
        c_ref, c_eff = complement(model_ref), complement(model_eff)
        if (c_ref, c_eff) == (pvar_ref, pvar_alt):
            return ("complement", pvar_alt)
        if (c_ref, c_eff) == (pvar_alt, pvar_ref):
            return ("complement_flipped", pvar_ref)
    return None


def load_model_targets(model_dir: Path, gene: str):
    targets = defaultdict(lambda: {"specs": set(), "tissues": set(), "genes": set(), "dbs": set()})
    for db_path in sorted(model_dir.glob("JTI_*.db")):
        tissue = db_path.stem.removeprefix("JTI_")
        with sqlite3.connect(db_path) as db:
            rows = db.execute(
                "select rsid,gene,ref_allele,eff_allele from weights "
                "where gene=? or gene like ?",
                (gene, gene + ".%"),
            ).fetchall()
        for rsid, model_gene, model_ref, model_eff in rows:
            target = targets[str(rsid)]
            target["specs"].add((str(model_ref).upper(), str(model_eff).upper()))
            target["tissues"].add(tissue)
            target["genes"].add(str(model_gene))
            target["dbs"].add(db_path)
    if not targets:
        raise RuntimeError(f"gene absent from JTI models: {gene}")
    return targets


def open_pvar(path: Path):
    if path.suffix == ".zst":
        process = subprocess.Popen(["zstdcat", str(path)], stdout=subprocess.PIPE, text=True)
        if process.stdout is None:
            raise RuntimeError("zstdcat did not provide stdout")
        return process, process.stdout
    return None, path.open()


def scan_pvar(path: Path, chromosome: str, wanted: set[str]):
    records = defaultdict(list)
    process, handle = open_pvar(path)
    try:
        header = None
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = [x.lstrip("#") for x in fields]
                continue
            row = dict(zip(header, fields))
            rsid = row.get("ID", "")
            if rsid not in wanted:
                continue
            chrom = row["CHROM"].removeprefix("chr")
            if chrom != str(chromosome).removeprefix("chr"):
                continue
            records[rsid].append(
                {
                    "chromosome": chrom,
                    "position": int(row["POS"]),
                    "ref": row["REF"].upper(),
                    "alts": tuple(x.upper() for x in row["ALT"].split(",")),
                }
            )
    finally:
        handle.close()
        if process is not None and process.wait() != 0:
            raise RuntimeError("zstdcat failed while reading PVAR")
    return records


def resolve_targets(model_targets, pvar_records):
    resolved, report = {}, []
    for rsid in sorted(model_targets):
        model = model_targets[rsid]
        base = {
            "gene": ",".join(sorted(model["genes"])),
            "rsid": rsid,
            "tissues": ",".join(sorted(model["tissues"])),
            "model_alleles": ",".join(f"{a}/{b}" for a, b in sorted(model["specs"])),
            "chromosome": "NA",
            "position": "NA",
            "pvar_ref": "NA",
            "pvar_alt": "NA",
            "selected_alt_index": "NA",
            "reference_effect_allele": "NA",
            "match_status": "NA",
            "multiallelic_reference": "NA",
            "missingness": "NA",
            "effect_allele_frequency": "NA",
            "variance": "NA",
            "final_status": "NA",
        }
        if not RSID.fullmatch(rsid):
            base["final_status"] = "missing_rsid"
            report.append(base)
            continue
        if any(is_palindromic(a, b) for a, b in model["specs"]):
            base["final_status"] = "palindrome_removed"
            report.append(base)
            continue
        per_spec = []
        for model_ref, model_eff in sorted(model["specs"]):
            candidates = []
            for record in pvar_records.get(rsid, []):
                for alt_index, alt in enumerate(record["alts"], 1):
                    matched = match_model_alleles(model_ref, model_eff, record["ref"], alt)
                    if matched and matched[0] != "palindrome_removed":
                        candidates.append((record, alt_index, alt, matched[0], matched[1]))
            if len(candidates) != 1:
                per_spec = []
                base["final_status"] = "missing_reference" if not candidates else "ambiguous"
                break
            per_spec.append(candidates[0])
        if not per_spec:
            if base["final_status"] == "NA":
                base["final_status"] = "allele_conflict"
            report.append(base)
            continue
        keys = {
            (
                x[0]["chromosome"],
                x[0]["position"],
                x[0]["ref"],
                x[2],
                x[1],
                x[4],
            )
            for x in per_spec
        }
        if len(keys) != 1:
            base["final_status"] = "model_effect_orientation_conflict"
            report.append(base)
            continue
        chromosome, position, ref, alt, alt_index, effect_allele = next(iter(keys))
        statuses = sorted({x[3] for x in per_spec})
        selected_record = per_spec[0][0]
        site_alts = {
            alt_value
            for record in pvar_records.get(rsid, [])
            if (
                record["chromosome"] == selected_record["chromosome"]
                and record["position"] == selected_record["position"]
                and record["ref"] == selected_record["ref"]
            )
            for alt_value in record["alts"]
        }
        alt_count = len(site_alts)
        base.update(
            chromosome=chromosome,
            position=position,
            pvar_ref=ref,
            pvar_alt=alt,
            selected_alt_index=alt_index,
            reference_effect_allele=effect_allele,
            match_status=",".join(statuses),
            multiallelic_reference=str(alt_count > 1).lower(),
            final_status="matched_reference",
        )
        resolved[rsid] = base
        report.append(base)
    return resolved, report


def load_fai(path: Path):
    index = {}
    with path.open() as handle:
        for line in handle:
            name, length, offset, line_bases, line_width, *_ = line.rstrip("\n").split("\t")
            index[name] = tuple(map(int, (length, offset, line_bases, line_width)))
    return index


def fasta_sequence(fasta: Path, fai, chromosome: str, position: int, length: int):
    contig = chromosome if chromosome in fai else "chr" + chromosome
    if contig not in fai:
        raise RuntimeError(f"FASTA contig absent: {chromosome}")
    contig_length, offset, line_bases, line_width = fai[contig]
    if position < 1 or position + length - 1 > contig_length:
        raise RuntimeError("FASTA query outside contig")
    result = bytearray()
    with fasta.open("rb") as handle:
        zero = position - 1
        handle.seek(offset + (zero // line_bases) * line_width + zero % line_bases)
        while len(result) < length:
            byte = handle.read(1)
            if not byte:
                raise RuntimeError("unexpected FASTA EOF")
            if byte not in (b"\n", b"\r"):
                result.extend(byte)
    return result.decode().upper()


def validate_reference(resolved, fasta: Path):
    fai_path = Path(str(fasta) + ".fai")
    if not fai_path.is_file():
        raise RuntimeError(f"FASTA index absent: {fai_path}")
    fai = load_fai(fai_path)
    bad = []
    for rsid, target in resolved.items():
        observed = fasta_sequence(
            fasta,
            fai,
            str(target["chromosome"]),
            int(target["position"]),
            len(str(target["pvar_ref"])),
        )
        if observed != target["pvar_ref"]:
            bad.append((rsid, target["pvar_ref"], observed))
    if bad:
        raise RuntimeError(f"reference_mismatch: {bad[:3]}")
    return fai_path


def run_plink_export(plink2: Path, pfile: Path, ids: Path, output_prefix: Path):
    command = [
        str(plink2),
        "--pfile",
        str(pfile),
        "--extract",
        str(ids),
        "--export",
        "vcf",
        "--out",
        str(output_prefix),
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError("PLINK2 VCF export failed: " + completed.stderr[-2000:])
    vcf = Path(str(output_prefix) + ".vcf")
    if not vcf.is_file():
        raise RuntimeError("PLINK2 did not create the expected VCF")
    return command, completed, vcf


def parse_effect_dosages(vcf: Path, resolved):
    dosages = {}
    with vcf.open() as handle:
        samples = []
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "#CHROM":
                samples = fields[9:]
                continue
            chrom, pos, rsid, ref, alt_field = fields[:5]
            if rsid not in resolved:
                continue
            target = resolved[rsid]
            alts = alt_field.upper().split(",")
            alt = str(target["pvar_alt"])
            if (
                chrom.removeprefix("chr") != str(target["chromosome"]).removeprefix("chr")
                or int(pos) != int(target["position"])
                or ref.upper() != target["pvar_ref"]
                or alt not in alts
            ):
                continue
            if rsid in dosages:
                raise RuntimeError(f"ambiguous exported VCF record: {rsid}")
            selected_alt_index = alts.index(alt) + 1
            effect_index = 0 if target["reference_effect_allele"] == ref.upper() else selected_alt_index
            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise RuntimeError("VCF GT field absent")
            gt_index = format_fields.index("GT")
            values = []
            for sample in fields[9:]:
                parts = sample.split(":")
                gt = parts[gt_index] if gt_index < len(parts) else "."
                alleles = re.split(r"[|/]", gt)
                if not alleles or any(x == "." for x in alleles):
                    values.append(float("nan"))
                else:
                    values.append(float(sum(int(x) == effect_index for x in alleles)))
            dosages[rsid] = np.asarray(values, dtype=float)
    missing = sorted(set(resolved) - set(dosages))
    if missing:
        raise RuntimeError(f"matched reference variants absent from exported VCF: {missing[:5]}")
    return samples, dosages


def qc_dosages(resolved, report, dosages, max_missing: float, min_maf: float):
    rows = {row["rsid"]: row for row in report}
    retained = {}
    for rsid in sorted(resolved):
        values = dosages[rsid]
        missingness = float(np.isnan(values).mean())
        observed = values[~np.isnan(values)]
        row = rows[rsid]
        row["missingness"] = missingness
        if observed.size == 0:
            row["final_status"] = "all_genotypes_missing"
            continue
        frequency = float(observed.mean() / 2.0)
        maf = min(frequency, 1.0 - frequency)
        row["effect_allele_frequency"] = frequency
        if missingness > max_missing:
            row["final_status"] = "missingness_removed"
            continue
        if maf < min_maf:
            row["final_status"] = "maf_removed"
            continue
        filled = values.copy()
        filled[np.isnan(filled)] = observed.mean()
        variance = float(np.var(filled, ddof=1))
        row["variance"] = variance
        if not math.isfinite(variance) or variance <= 0:
            row["final_status"] = "zero_variance_removed"
            continue
        row["final_status"] = "retained"
        retained[rsid] = filled
    if not retained:
        raise RuntimeError("no model variants survived reference matching and genotype QC")
    return retained


def covariance_matrix(retained):
    rsids = sorted(retained)
    matrix = np.atleast_2d(np.cov(np.vstack([retained[x] for x in rsids]), ddof=1))
    if matrix.shape != (len(rsids), len(rsids)) or not np.isfinite(matrix).all():
        raise RuntimeError("nonfinite or malformed covariance matrix")
    if not np.allclose(matrix, matrix.T, atol=1e-12, rtol=1e-10):
        raise RuntimeError("covariance matrix is asymmetric")
    eigenvalues = np.linalg.eigvalsh(matrix)
    tolerance = 1e-10 * max(1.0, float(np.max(eigenvalues)))
    if float(np.min(eigenvalues)) < -tolerance:
        raise RuntimeError("covariance matrix is materially indefinite")
    return rsids, matrix, eigenvalues


def write_covariance(path: Path, gene: str, rsids, matrix):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["GENE", "RSID1", "RSID2", "VALUE"])
                for i, rsid1 in enumerate(rsids):
                    for j, rsid2 in enumerate(rsids):
                        writer.writerow([gene, rsid1, rsid2, format(float(matrix[i, j]), ".17g")])


def write_tsv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["rsid"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def normalized_gwas_rsids(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "rsid" not in (reader.fieldnames or []):
            raise RuntimeError("normalized GWAS lacks rsid column")
        return {row["rsid"] for row in reader if RSID.fullmatch(row.get("rsid", ""))}


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--models-folder", type=Path, required=True)
    result.add_argument("--gene", required=True, help="unversioned or exact Ensembl gene ID")
    result.add_argument("--chromosome", required=True)
    result.add_argument("--pfile", type=Path, required=True, help="PLINK2 prefix without extension")
    result.add_argument("--plink2", type=Path, required=True)
    result.add_argument("--fasta", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--matching-report", type=Path, required=True)
    result.add_argument("--method-json", type=Path, required=True)
    result.add_argument("--work-dir", type=Path, required=True)
    result.add_argument("--normalized-gwas", type=Path)
    result.add_argument("--cleared-output", type=Path)
    result.add_argument("--max-missing", type=float, default=0.02)
    result.add_argument("--min-maf", type=float, default=0.0)
    return result


def main():
    args = parser().parse_args()
    if args.work_dir.exists() and any(args.work_dir.iterdir()):
        raise RuntimeError(f"work directory must be new or empty: {args.work_dir}")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    pvar = Path(str(args.pfile) + ".pvar")
    if not pvar.is_file():
        pvar = Path(str(args.pfile) + ".pvar.zst")
    for required in (pvar, Path(str(args.pfile) + ".pgen"), Path(str(args.pfile) + ".psam"), args.fasta, args.plink2):
        if not required.exists():
            raise RuntimeError(f"required input absent: {required}")

    model_targets = load_model_targets(args.models_folder, args.gene)
    pvar_records = scan_pvar(pvar, args.chromosome, set(model_targets))
    resolved, report = resolve_targets(model_targets, pvar_records)
    fai = validate_reference(resolved, args.fasta)
    ids = args.work_dir / "extract-rsids.txt"
    ids.write_text("".join(f"{x}\n" for x in sorted(resolved)))
    command, completed, vcf = run_plink_export(args.plink2, args.pfile, ids, args.work_dir / "reference")
    samples, dosages = parse_effect_dosages(vcf, resolved)
    retained = qc_dosages(resolved, report, dosages, args.max_missing, args.min_maf)
    rsids, matrix, eigenvalues = covariance_matrix(retained)
    write_covariance(args.output, args.gene, rsids, matrix)
    write_tsv(args.matching_report, report)

    cleared_count = None
    if bool(args.normalized_gwas) != bool(args.cleared_output):
        raise RuntimeError("--normalized-gwas and --cleared-output must be supplied together")
    if args.normalized_gwas:
        intersection = sorted(set(rsids) & normalized_gwas_rsids(args.normalized_gwas))
        write_tsv(args.cleared_output, [{"rsid": x} for x in intersection])
        cleared_count = len(intersection)

    contributing_dbs = sorted({p for target in model_targets.values() for p in target["dbs"]})
    method = {
        "method": "metaxcan-common-reference-dosage-covariance-v1",
        "gene": args.gene,
        "chromosome": str(args.chromosome),
        "reference_build": "GRCh37/hg19",
        "multiallelic_mode": "exact_alt_match",
        "indel_normalization": "assumed_pre_normalized_and_exact_string_matched",
        "palindromic_policy": "remove",
        "strand_policy": "allow_nonpalindromic_snp_complement",
        "dosage_orientation": "model_effect_allele",
        "missing_genotype_policy": "mean_impute_after_missingness_qc",
        "covariance_ddof": 1,
        "max_missing": args.max_missing,
        "min_maf": args.min_maf,
        "n_model_rsids": len(model_targets),
        "n_reference_matched": len(resolved),
        "n_retained": len(rsids),
        "n_samples": len(samples),
        "n_cleared_for_gwas": cleared_count,
        "min_eigenvalue": float(np.min(eigenvalues)),
        "max_eigenvalue": float(np.max(eigenvalues)),
        "plink_command": command,
        "plink_stdout_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
        "plink_stderr_sha256": hashlib.sha256(completed.stderr.encode()).hexdigest(),
        "inputs": {
            "pgen": str(Path(str(args.pfile) + ".pgen").resolve()),
            "pvar": str(pvar.resolve()),
            "psam": str(Path(str(args.pfile) + ".psam").resolve()),
            "fasta": str(args.fasta.resolve()),
            "fasta_fai_sha256": sha256(fai),
            "model_databases": [{"path": str(p.resolve()), "sha256": sha256(p)} for p in contributing_dbs],
            "normalized_gwas": str(args.normalized_gwas.resolve()) if args.normalized_gwas else None,
        },
        "outputs": {
            "covariance": {"path": str(args.output.resolve()), "sha256": sha256(args.output)},
            "matching_report": {"path": str(args.matching_report.resolve()), "sha256": sha256(args.matching_report)},
            "cleared_snps": {"path": str(args.cleared_output.resolve()), "sha256": sha256(args.cleared_output)} if args.cleared_output else None,
        },
    }
    args.method_json.parent.mkdir(parents=True, exist_ok=True)
    args.method_json.write_text(json.dumps(method, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
