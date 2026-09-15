#!/usr/bin/env python3
"""Build an audited genome-wide JTI S-MultiXcan dosage-covariance resource.

The implementation indexes all JTI models once, scans each hg19 PVAR once,
exports one genotype VCF per chromosome, and reuses each dosage vector across
all genes.  It deliberately does not merge tissue covariance files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import itertools
import json
import math
import os
import re
import sqlite3
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from smultixcan_genome_core import normalized_gene, sha256  # noqa: E402


REPORT_FIELDS = [
    "gene", "model_gene_entities", "rsid", "tissues", "model_alleles",
    "chromosome", "position", "pvar_ref", "pvar_alt",
    "selected_alt_index", "reference_effect_allele", "match_status",
    "multiallelic_reference", "missingness", "effect_allele_frequency",
    "variance", "final_status",
]
GENE_QC_FIELDS = [
    "gene", "chromosome", "build_status", "n_model_rsids",
    "n_reference_matched", "n_retained", "n_samples", "min_eigenvalue",
    "max_eigenvalue", "psd_tolerance", "covariance_cells",
]
RSID = re.compile(r"^rs[0-9]+$")


def load_matching_module(project_root: Path):
    del project_root  # Kept in the CLI for project provenance and working-directory validation.
    path = HERE / "build_smultixcan_covariance.py"
    if not path.is_file():
        raise RuntimeError(f"single-gene covariance implementation absent: {path}")
    spec = importlib.util.spec_from_file_location("m1_single_gene_covariance", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module, path


def deterministic_gzip_text(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(compressed, newline="")
    return raw, compressed, text


def parse_gene_filter(values):
    result = set()
    for value in values or []:
        result.update(normalized_gene(x.strip()) for x in value.split(",") if x.strip())
    return result


def initialize_index(path: Path, model_dir: Path, gene_filter) -> dict:
    if path.exists():
        path.unlink()
    connection = sqlite3.connect(path)
    connection.execute("pragma journal_mode=wal")
    connection.execute("pragma synchronous=normal")
    connection.execute("pragma temp_store=memory")
    connection.executescript(
        """
        create table model(
          gene text not null, model_gene text not null, rsid text not null,
          ref text not null, eff text not null, tissue text not null,
          primary key(gene,rsid,ref,eff,tissue,model_gene)
        ) without rowid;
        create table gene_name(
          gene text not null, gene_name text not null,
          primary key(gene,gene_name)
        ) without rowid;
        create table gene_state(
          gene text primary key, chromosome text, match_status text not null,
          n_model_rsids integer not null, n_reference_matched integer not null
        ) without rowid;
        create table resolved(
          gene text not null, rsid text not null, chromosome text not null,
          position integer not null, pvar_ref text not null, pvar_alt text not null,
          selected_alt_index integer not null, effect_allele text not null,
          primary key(gene,rsid)
        ) without rowid;
        create table report(
          gene text not null, model_gene_entities text not null, rsid text not null,
          tissues text not null, model_alleles text not null, chromosome text not null,
          position text not null, pvar_ref text not null, pvar_alt text not null,
          selected_alt_index text not null, reference_effect_allele text not null,
          match_status text not null, multiallelic_reference text not null,
          missingness text not null, effect_allele_frequency text not null,
          variance text not null, final_status text not null,
          primary key(gene,rsid)
        ) without rowid;
        """
    )
    db_meta = []
    model_gene_entities = defaultdict(set)
    n_weight_rows = 0
    for db_path in sorted(model_dir.glob("JTI_*.db")):
        tissue = db_path.stem.removeprefix("JTI_")
        with sqlite3.connect(db_path) as source:
            weight_rows = source.execute(
                "select rsid,gene,ref_allele,eff_allele from weights"
            )
            batch = []
            for rsid, model_gene, ref, eff in weight_rows:
                gene = normalized_gene(model_gene)
                if gene_filter and gene not in gene_filter:
                    continue
                model_gene_entities[gene].add(str(model_gene))
                batch.append((gene, str(model_gene), str(rsid), str(ref).upper(), str(eff).upper(), tissue))
                if len(batch) >= 50_000:
                    connection.executemany("insert or ignore into model values(?,?,?,?,?,?)", batch)
                    n_weight_rows += len(batch)
                    batch.clear()
            if batch:
                connection.executemany("insert or ignore into model values(?,?,?,?,?,?)", batch)
                n_weight_rows += len(batch)
            for model_gene, gene_name in source.execute("select gene,genename from extra"):
                gene = normalized_gene(model_gene)
                if gene_filter and gene not in gene_filter:
                    continue
                connection.execute(
                    "insert or ignore into gene_name values(?,?)",
                    (gene, str(gene_name) if gene_name is not None else "NA"),
                )
        connection.commit()
        db_meta.append({"path": str(db_path.resolve()), "sha256": sha256(db_path), "tissue": tissue})
    collision = {gene: sorted(values) for gene, values in model_gene_entities.items() if len(values) > 1}
    if collision:
        raise RuntimeError(f"Ensembl version trimming collision: {list(collision.items())[:3]}")
    connection.execute("create index model_rsid on model(rsid)")
    connection.commit()
    counts = {
        "model_databases": len(db_meta),
        "indexed_weight_rows": n_weight_rows,
        "genes": connection.execute("select count(distinct gene) from model").fetchone()[0],
        "unique_gene_rsids": connection.execute("select count(*) from (select distinct gene,rsid from model)").fetchone()[0],
        "unique_rsids": connection.execute("select count(distinct rsid) from model").fetchone()[0],
    }
    connection.close()
    return {"model_databases": db_meta, "counts": counts}


def scan_pvars(index: Path, ld_dir: Path, chromosomes, matching):
    with sqlite3.connect(index) as connection:
        wanted = {row[0] for row in connection.execute("select distinct rsid from model")}
    records = defaultdict(list)
    pvar_meta = []
    for chromosome in chromosomes:
        pvar = ld_dir / f"eur_chr{chromosome}.pvar"
        if not pvar.is_file():
            pvar = Path(str(pvar) + ".zst")
        if not pvar.is_file():
            raise RuntimeError(f"PVAR absent: {pvar}")
        process, handle = matching.open_pvar(pvar)
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
                if chrom != str(chromosome):
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
                raise RuntimeError(f"zstdcat failed: {pvar}")
        pvar_meta.append({"chromosome": str(chromosome), "path": str(pvar.resolve()), "sha256": sha256(pvar)})
    return records, pvar_meta


class ReferenceReader:
    def __init__(self, fasta: Path, matching):
        self.fasta = fasta
        self.fai_path = Path(str(fasta) + ".fai")
        if not self.fai_path.is_file():
            raise RuntimeError(f"FASTA index absent: {self.fai_path}")
        self.index = matching.load_fai(self.fai_path)
        self.handle = fasta.open("rb")

    def sequence(self, chromosome, position, length):
        contig = str(chromosome) if str(chromosome) in self.index else "chr" + str(chromosome)
        if contig not in self.index:
            raise RuntimeError(f"FASTA contig absent: {chromosome}")
        contig_length, offset, line_bases, line_width = self.index[contig]
        if position < 1 or position + length - 1 > contig_length:
            raise RuntimeError("FASTA query outside contig")
        zero = position - 1
        self.handle.seek(offset + (zero // line_bases) * line_width + zero % line_bases)
        result = bytearray()
        while len(result) < length:
            byte = self.handle.read(1)
            if not byte:
                raise RuntimeError("unexpected FASTA EOF")
            if byte not in (b"\n", b"\r"):
                result.extend(byte)
        return result.decode().upper()

    def close(self):
        self.handle.close()


def grouped_model_targets(connection):
    cursor = connection.execute(
        "select gene,model_gene,rsid,ref,eff,tissue from model order by gene,rsid,ref,eff,tissue"
    )
    for gene, rows in itertools.groupby(cursor, key=lambda row: row[0]):
        targets = defaultdict(lambda: {"specs": set(), "tissues": set(), "genes": set(), "dbs": set()})
        for _, model_gene, rsid, ref, eff, tissue in rows:
            target = targets[rsid]
            target["specs"].add((ref, eff))
            target["tissues"].add(tissue)
            target["genes"].add(model_gene)
        yield gene, targets


def build_matching_index(index: Path, pvar_records, fasta: Path, matching):
    connection = sqlite3.connect(index)
    reader = ReferenceReader(fasta, matching)
    ref_cache = {}
    try:
        for number, (gene, targets) in enumerate(grouped_model_targets(connection), 1):
            resolved, report = matching.resolve_targets(targets, {r: pvar_records.get(r, []) for r in targets})
            bad_reference = set()
            for rsid, target in resolved.items():
                key = (str(target["chromosome"]), int(target["position"]), str(target["pvar_ref"]))
                if key not in ref_cache:
                    ref_cache[key] = reader.sequence(key[0], key[1], len(key[2])) == key[2]
                if not ref_cache[key]:
                    bad_reference.add(rsid)
            for row in report:
                if row["rsid"] in bad_reference:
                    row["final_status"] = "reference_mismatch"
                    resolved.pop(row["rsid"], None)
            chroms = {str(x["chromosome"]) for x in resolved.values()}
            if len(chroms) > 1:
                state = "gene_cross_chromosome"
                chromosome = "NA"
                for row in report:
                    if row["rsid"] in resolved:
                        row["final_status"] = state
                resolved = {}
            elif len(chroms) == 1:
                state = "matched_reference"
                chromosome = next(iter(chroms))
            else:
                state = "no_reference_matched_variants"
                chromosome = "NA"
            connection.execute(
                "insert into gene_state values(?,?,?,?,?)",
                (gene, chromosome, state, len(targets), len(resolved)),
            )
            for rsid, target in resolved.items():
                connection.execute(
                    "insert into resolved values(?,?,?,?,?,?,?,?)",
                    (
                        gene, rsid, str(target["chromosome"]), int(target["position"]),
                        str(target["pvar_ref"]), str(target["pvar_alt"]),
                        int(target["selected_alt_index"]), str(target["reference_effect_allele"]),
                    ),
                )
            report_rows = []
            for row in report:
                report_rows.append(
                    (
                        gene, row["gene"], row["rsid"], row["tissues"], row["model_alleles"],
                        str(row["chromosome"]), str(row["position"]), str(row["pvar_ref"]),
                        str(row["pvar_alt"]), str(row["selected_alt_index"]),
                        str(row["reference_effect_allele"]), str(row["match_status"]),
                        str(row["multiallelic_reference"]), str(row["missingness"]),
                        str(row["effect_allele_frequency"]), str(row["variance"]),
                        str(row["final_status"]),
                    )
                )
            connection.executemany("insert into report values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", report_rows)
            if number % 250 == 0:
                connection.commit()
        connection.commit()
    finally:
        reader.close()
        connection.close()


def write_extract(path: Path, rsids):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{x}\n" for x in sorted(rsids)))


def run_plink(plink2: Path, pfile: Path, extract: Path, prefix: Path):
    command = [
        str(plink2), "--pfile", str(pfile), "--extract", str(extract),
        "--export", "vcf", "--out", str(prefix),
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    (prefix.parent / "plink.stdout.txt").write_text(completed.stdout)
    (prefix.parent / "plink.stderr.txt").write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError(f"PLINK2 chromosome export failed: {completed.stderr[-2000:]}")
    vcf = Path(str(prefix) + ".vcf")
    if not vcf.is_file():
        raise RuntimeError(f"PLINK2 VCF absent: {vcf}")
    return command, completed, vcf


def target_key(row):
    return (
        row["rsid"], str(row["chromosome"]), int(row["position"]),
        row["pvar_ref"], row["pvar_alt"], row["effect_allele"],
    )


def parse_batch_dosages(vcf: Path, rows):
    targets = defaultdict(list)
    for row in rows:
        key = target_key(row)
        if key not in targets[row["rsid"]]:
            targets[row["rsid"]].append(key)
    dosages = {}
    samples = []
    with vcf.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "#CHROM":
                samples = fields[9:]
                continue
            chromosome, position, rsid, ref, alt_field = fields[:5]
            if rsid not in targets:
                continue
            alts = alt_field.upper().split(",")
            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise RuntimeError("VCF GT field absent")
            gt_index = format_fields.index("GT")
            parsed = []
            for sample in fields[9:]:
                parts = sample.split(":")
                gt = parts[gt_index] if gt_index < len(parts) else "."
                alleles = re.split(r"[|/]", gt)
                parsed.append(None if not alleles or any(x == "." for x in alleles) else tuple(map(int, alleles)))
            for key in targets[rsid]:
                _, target_chrom, target_pos, target_ref, target_alt, effect = key
                if (
                    chromosome.removeprefix("chr") != target_chrom.removeprefix("chr")
                    or int(position) != target_pos or ref.upper() != target_ref
                    or target_alt not in alts
                ):
                    continue
                if key in dosages:
                    raise RuntimeError(f"ambiguous exported VCF target: {key}")
                effect_index = 0 if effect == target_ref else alts.index(target_alt) + 1
                dosages[key] = np.asarray(
                    [float("nan") if gt is None else float(sum(x == effect_index for x in gt)) for gt in parsed],
                    dtype=float,
                )
    return samples, dosages


def qc_dosage(values, max_missing, min_maf):
    if values is None:
        return None, {"status": "vcf_absent", "missingness": "NA", "frequency": "NA", "variance": "NA"}
    missingness = float(np.isnan(values).mean())
    observed = values[~np.isnan(values)]
    if observed.size == 0:
        return None, {"status": "all_genotypes_missing", "missingness": missingness, "frequency": "NA", "variance": "NA"}
    frequency = float(observed.mean() / 2.0)
    maf = min(frequency, 1.0 - frequency)
    if missingness > max_missing:
        return None, {"status": "missingness_removed", "missingness": missingness, "frequency": frequency, "variance": "NA"}
    if maf < min_maf:
        return None, {"status": "maf_removed", "missingness": missingness, "frequency": frequency, "variance": "NA"}
    filled = values.copy()
    filled[np.isnan(filled)] = observed.mean()
    variance = float(np.var(filled, ddof=1))
    if not math.isfinite(variance) or variance <= 0:
        return None, {"status": "zero_variance_removed", "missingness": missingness, "frequency": frequency, "variance": variance}
    return filled, {"status": "retained", "missingness": missingness, "frequency": frequency, "variance": variance}


def chromosome_rows(connection, chromosome):
    connection.row_factory = sqlite3.Row
    return list(
        connection.execute(
            "select r.gene,r.rsid,r.chromosome,r.position,r.pvar_ref,r.pvar_alt,r.selected_alt_index,r.effect_allele,"
            "g.n_model_rsids,g.n_reference_matched "
            "from resolved r join gene_state g on g.gene=r.gene "
            "where r.chromosome=? order by r.gene,r.rsid",
            (str(chromosome),),
        )
    )


def build_chromosome_shard(args, chromosome, index, matching, include_header):
    shard_dir = args.output_dir / "shards" / f"chr{chromosome}"
    work_dir = args.output_dir / "work" / f"chr{chromosome}"
    shard_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    marker = shard_dir / "method.json"
    if args.resume and marker.is_file():
        method = json.loads(marker.read_text())
        outputs = method.get("outputs", {})
        if all(Path(v["path"]).is_file() and sha256(Path(v["path"])) == v["sha256"] for v in outputs.values()):
            return method
    with sqlite3.connect(index) as connection:
        rows = chromosome_rows(connection, chromosome)
    if not rows:
        raise RuntimeError(f"no resolved model variants on chromosome {chromosome}")
    rsids = {row["rsid"] for row in rows}
    extract = work_dir / "extract-rsids.txt"
    write_extract(extract, rsids)
    pfile = args.ld_reference_folder / f"eur_chr{chromosome}"
    command, completed, vcf = run_plink(args.plink2, pfile, extract, work_dir / "reference")
    samples, dosage_raw = parse_batch_dosages(vcf, rows)
    key_qc = {}
    retained = {}
    for row in rows:
        key = target_key(row)
        if key in key_qc:
            continue
        values, qc = qc_dosage(dosage_raw.get(key), args.max_missing, args.min_maf)
        key_qc[key] = qc
        if values is not None:
            retained[key] = values
    with sqlite3.connect(index) as connection:
        for row in rows:
            qc = key_qc[target_key(row)]
            connection.execute(
                "update report set missingness=?,effect_allele_frequency=?,variance=?,final_status=? where gene=? and rsid=?",
                (
                    str(qc["missingness"]), str(qc["frequency"]), str(qc["variance"]),
                    qc["status"], row["gene"], row["rsid"],
                ),
            )
        connection.commit()
    covariance = shard_dir / "covariance.txt.gz"
    covariance_tmp = shard_dir / "covariance.txt.gz.partial"
    gene_qc = shard_dir / "gene_qc.tsv.gz"
    retained_file = shard_dir / "retained_rsids.tsv.gz"
    raw, compressed, text = deterministic_gzip_text(covariance_tmp)
    writer = csv.writer(text, delimiter="\t", lineterminator="\n")
    if include_header:
        writer.writerow(["GENE", "RSID1", "RSID2", "VALUE"])
    qc_rows = []
    retained_rsids = set()
    for gene, gene_rows_iter in itertools.groupby(rows, key=lambda row: row["gene"]):
        gene_rows = list(gene_rows_iter)
        vectors = {row["rsid"]: retained[target_key(row)] for row in gene_rows if target_key(row) in retained}
        n_model = int(gene_rows[0]["n_model_rsids"])
        n_reference = int(gene_rows[0]["n_reference_matched"])
        if not vectors:
            qc_rows.append({
                "gene": gene, "chromosome": str(chromosome), "build_status": "no_variants_after_genotype_qc",
                "n_model_rsids": n_model, "n_reference_matched": n_reference, "n_retained": 0,
                "n_samples": len(samples), "min_eigenvalue": "NA", "max_eigenvalue": "NA",
                "psd_tolerance": "NA", "covariance_cells": 0,
            })
            continue
        gene_rsids, matrix, eigenvalues = matching.covariance_matrix(vectors)
        max_eigen = float(np.max(eigenvalues))
        tolerance = 1e-10 * max(1.0, max_eigen)
        for i, rsid1 in enumerate(gene_rsids):
            for j, rsid2 in enumerate(gene_rsids):
                writer.writerow([gene, rsid1, rsid2, format(float(matrix[i, j]), ".17g")])
        retained_rsids.update(gene_rsids)
        qc_rows.append({
            "gene": gene, "chromosome": str(chromosome), "build_status": "success",
            "n_model_rsids": n_model, "n_reference_matched": n_reference,
            "n_retained": len(gene_rsids), "n_samples": len(samples),
            "min_eigenvalue": format(float(np.min(eigenvalues)), ".17g"),
            "max_eigenvalue": format(max_eigen, ".17g"),
            "psd_tolerance": format(tolerance, ".17g"),
            "covariance_cells": len(gene_rsids) ** 2,
        })
    text.flush()
    text.close()
    compressed.close()
    raw.close()
    os.replace(covariance_tmp, covariance)
    with gzip.open(gene_qc, "wt", newline="") as handle:
        out = csv.DictWriter(handle, fieldnames=GENE_QC_FIELDS, delimiter="\t", lineterminator="\n")
        out.writeheader(); out.writerows(qc_rows)
    with gzip.open(retained_file, "wt", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["rsid"]); out.writerows((x,) for x in sorted(retained_rsids))
    method = {
        "chromosome": str(chromosome), "status": "pass", "n_samples": len(samples),
        "n_resolved_rows": len(rows), "n_extract_rsids": len(rsids),
        "n_retained_rsids": len(retained_rsids), "n_genes": len(qc_rows),
        "n_successful_gene_covariances": sum(x["build_status"] == "success" for x in qc_rows),
        "plink_command": command, "plink_exit_code": completed.returncode,
        "outputs": {
            "covariance": {"path": str(covariance.resolve()), "sha256": sha256(covariance)},
            "gene_qc": {"path": str(gene_qc.resolve()), "sha256": sha256(gene_qc)},
            "retained_rsids": {"path": str(retained_file.resolve()), "sha256": sha256(retained_file)},
        },
    }
    marker.write_text(json.dumps(method, indent=2, sort_keys=True) + "\n")
    return method


def export_matching(index: Path, output: Path):
    with sqlite3.connect(index) as connection, gzip.open(output, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(REPORT_FIELDS)
        writer.writerows(
            connection.execute(
                "select gene,model_gene_entities,rsid,tissues,model_alleles,chromosome,position,pvar_ref,pvar_alt,"
                "selected_alt_index,reference_effect_allele,match_status,multiallelic_reference,missingness,"
                "effect_allele_frequency,variance,final_status from report order by gene,rsid"
            )
        )


def combine_outputs(args, chromosomes, shard_methods, index, index_meta, pvar_meta, implementation_path):
    covariance = args.output_dir / "snp_covariance.txt.gz"
    with covariance.open("wb") as output:
        for chromosome in chromosomes:
            source = args.output_dir / "shards" / f"chr{chromosome}" / "covariance.txt.gz"
            with source.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    output.write(chunk)
    retained = set()
    gene_qc_rows = []
    for chromosome in chromosomes:
        shard = args.output_dir / "shards" / f"chr{chromosome}"
        with gzip.open(shard / "retained_rsids.tsv.gz", "rt", newline="") as handle:
            retained.update(row["rsid"] for row in csv.DictReader(handle, delimiter="\t"))
        with gzip.open(shard / "gene_qc.tsv.gz", "rt", newline="") as handle:
            gene_qc_rows.extend(csv.DictReader(handle, delimiter="\t"))
    retained_path = args.output_dir / "retained_rsids.tsv.gz"
    with gzip.open(retained_path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["rsid"]); writer.writerows((x,) for x in sorted(retained))
    gene_qc_path = args.output_dir / "gene_qc.tsv.gz"
    with gzip.open(gene_qc_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GENE_QC_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(sorted(gene_qc_rows, key=lambda x: x["gene"]))
    matching_path = args.output_dir / "variant_matching.tsv.gz"
    export_matching(index, matching_path)
    with sqlite3.connect(index) as connection:
        state_counts = dict(connection.execute("select match_status,count(*) from gene_state group by match_status"))
        report_counts = dict(connection.execute("select final_status,count(*) from report group by final_status"))
        n_gene_states = connection.execute("select count(*) from gene_state").fetchone()[0]
        n_report_rows = connection.execute("select count(*) from report").fetchone()[0]
    successful = [x for x in gene_qc_rows if x["build_status"] == "success"]
    scope = "genome_wide" if not args.gene_filter and set(chromosomes) == {str(x) for x in range(1, 23)} else "test_subset"
    validation = {
        "status": "pass" if successful and all(x["status"] == "pass" for x in shard_methods) else "fail",
        "scope": scope,
        "checks": {
            "all_chromosome_shards_pass": all(x["status"] == "pass" for x in shard_methods),
            "covariance_exists_nonempty": covariance.is_file() and covariance.stat().st_size > 0,
            "retained_rsid_set_nonempty": bool(retained),
            "successful_gene_covariance_set_nonempty": bool(successful),
            "all_successful_gene_matrices_psd": all(
                float(x["min_eigenvalue"]) >= -float(x["psd_tolerance"]) for x in successful
            ),
            "all_reference_sample_counts_identical": len({x["n_samples"] for x in successful}) == 1,
            "model_gene_state_closure": n_gene_states == index_meta["counts"]["genes"],
            "model_gene_rsid_report_closure": n_report_rows == index_meta["counts"]["unique_gene_rsids"],
            "matched_gene_shard_closure": len(gene_qc_rows) == state_counts.get("matched_reference", 0),
        },
        "counts": {
            **index_meta["counts"], "retained_rsids": len(retained),
            "gene_covariance_rows": len(gene_qc_rows),
            "successful_gene_covariances": len(successful),
            "ordered_covariance_cells": sum(int(x["covariance_cells"]) for x in successful),
        },
        "gene_match_status_counts": state_counts,
        "variant_final_status_counts": report_counts,
    }
    validation_path = args.output_dir / "validation-report.json"
    validation_path.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    method = {
        "method": "genome-wide-common-reference-dosage-covariance-v1",
        "scope": scope,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "reference_build": "GRCh37/hg19", "ancestry": "1000_Genomes_EUR",
        "chromosomes": list(map(str, chromosomes)),
        "multiallelic_mode": "exact_alt_match", "palindromic_policy": "remove",
        "strand_policy": "allow_nonpalindromic_snp_complement",
        "dosage_orientation": "model_effect_allele", "covariance_ddof": 1,
        "missing_genotype_policy": "mean_impute_after_missingness_qc",
        "max_missing": args.max_missing, "min_maf": args.min_maf,
        "gene_filter": sorted(args.gene_filter) if args.gene_filter else None,
        "inputs": {
            "model_databases": index_meta["model_databases"], "pvars": pvar_meta,
            "fasta": str(args.fasta.resolve()), "fasta_sha256": sha256(args.fasta),
            "fasta_fai_sha256": sha256(Path(str(args.fasta) + ".fai")),
            "single_gene_matching_implementation": {
                "path": str(implementation_path.resolve()), "sha256": sha256(implementation_path),
            },
        },
        "chromosome_shards": shard_methods,
        "outputs": {
            "covariance": {"path": str(covariance.resolve()), "sha256": sha256(covariance)},
            "retained_rsids": {"path": str(retained_path.resolve()), "sha256": sha256(retained_path)},
            "gene_qc": {"path": str(gene_qc_path.resolve()), "sha256": sha256(gene_qc_path)},
            "variant_matching": {"path": str(matching_path.resolve()), "sha256": sha256(matching_path)},
            "validation_report": {"path": str(validation_path.resolve()), "sha256": sha256(validation_path)},
        },
    }
    (args.output_dir / "method.json").write_text(json.dumps(method, indent=2, sort_keys=True) + "\n")
    if validation["status"] != "pass":
        raise RuntimeError("genome-wide covariance validation failed")


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--project-root", type=Path, required=True)
    result.add_argument("--models-folder", type=Path, required=True)
    result.add_argument("--ld-reference-folder", type=Path, required=True)
    result.add_argument("--plink2", type=Path, required=True)
    result.add_argument("--fasta", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    result.add_argument("--chromosomes", nargs="+", default=[str(x) for x in range(1, 23)])
    result.add_argument("--genes", action="append", default=[], help="Optional comma-separated test subset")
    result.add_argument("--max-missing", type=float, default=0.02)
    result.add_argument("--min-maf", type=float, default=0.0)
    result.add_argument("--resume", action="store_true")
    return result


def main():
    args = parser().parse_args()
    args.project_root = args.project_root.resolve()
    args.models_folder = args.models_folder.resolve()
    args.ld_reference_folder = args.ld_reference_folder.resolve()
    args.plink2 = args.plink2.resolve()
    args.fasta = args.fasta.resolve()
    args.output_dir = args.output_dir.resolve()
    args.gene_filter = parse_gene_filter(args.genes)
    chromosomes = [str(x).removeprefix("chr") for x in args.chromosomes]
    if not chromosomes or any(not x.isdigit() or not 1 <= int(x) <= 22 for x in chromosomes):
        raise RuntimeError("chromosomes must be integers 1-22")
    if len(chromosomes) != len(set(chromosomes)):
        raise RuntimeError("duplicate chromosomes")
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.resume:
        raise RuntimeError(f"output directory must be new/empty unless --resume: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for path in (args.models_folder, args.ld_reference_folder, args.plink2, args.fasta, Path(str(args.fasta) + ".fai")):
        if not path.exists():
            raise RuntimeError(f"required input absent: {path}")
    matching, implementation_path = load_matching_module(args.project_root)
    current_parameters = {
        "chromosomes": chromosomes,
        "gene_filter": sorted(args.gene_filter),
        "max_missing": args.max_missing,
        "min_maf": args.min_maf,
        "batch_builder_sha256": sha256(Path(__file__)),
        "single_gene_implementation_sha256": sha256(implementation_path),
    }
    index = args.output_dir / "work" / "model-index.sqlite"
    index.parent.mkdir(parents=True, exist_ok=True)
    index_meta_path = args.output_dir / "work" / "index-method.json"
    if args.resume and index.is_file() and index_meta_path.is_file():
        index_meta = json.loads(index_meta_path.read_text())
        if index_meta.get("parameters") != current_parameters:
            raise RuntimeError("--resume parameters or implementation differ from the indexed run")
        for item in index_meta["model_databases"] + index_meta["pvars"]:
            path = Path(item["path"])
            if not path.is_file() or sha256(path) != item["sha256"]:
                raise RuntimeError(f"--resume input hash mismatch: {path}")
        if index_meta.get("fasta_sha256") != sha256(args.fasta) or index_meta.get("fasta_fai_sha256") != sha256(Path(str(args.fasta) + ".fai")):
            raise RuntimeError("--resume FASTA identity mismatch")
    else:
        index_meta = initialize_index(index, args.models_folder, args.gene_filter)
        pvar_records, pvar_meta = scan_pvars(index, args.ld_reference_folder, chromosomes, matching)
        build_matching_index(index, pvar_records, args.fasta, matching)
        index_meta["pvars"] = pvar_meta
        index_meta["parameters"] = current_parameters
        index_meta["fasta_sha256"] = sha256(args.fasta)
        index_meta["fasta_fai_sha256"] = sha256(Path(str(args.fasta) + ".fai"))
        index_meta_path.write_text(json.dumps(index_meta, indent=2, sort_keys=True) + "\n")
    pvar_meta = index_meta["pvars"]
    shard_methods = []
    for i, chromosome in enumerate(chromosomes):
        shard_methods.append(build_chromosome_shard(args, chromosome, index, matching, include_header=(i == 0)))
    combine_outputs(args, chromosomes, shard_methods, index, index_meta, pvar_meta, implementation_path)
    print(json.dumps({"status": "pass", "output_dir": str(args.output_dir)}, sort_keys=True))


if __name__ == "__main__":
    main()
