#!/usr/bin/env python3
"""Independent numeric verification of the frozen Module 3 v6 exploratory results.

Written 2026-09-10 for pre-submission audit B-12. validate_v5.py checks the
structural contract of the v5/v6 runs; this script independently re-derives the
key quantities promised by RESEARCH_PLAN §7.6 without trusting any stored
summary beyond the frozen inputs themselves:

1. input_hashes           every input_sources / locked_outputs sha256 in the
                          registry matches the file on disk.
2. universe_contract      8,101 rows, unique HGNC ids and symbols, one-to-one
                          Ensembl mapping across the PE(JTI) and GSE103927 sides.
3. strict_sets_rebuilt    the four PRIMARY_STRICT sets are rebuilt from the RAW
                          frozen sources (module-2 v2 ALT statistics; frozen PE
                          Whole Blood TWAS table), not from the result tables.
4. observed_intersections the four quadrant overlaps are rebuilt from those sets
                          and match directional_overlap.tsv; the only nonzero
                          primary intersection is WARS1 (alt_up_pe_negative).
5. fisher_recomputed      one-sided Fisher P is recomputed from the
                          hypergeometric tail for all 12 scenario-quadrant rows.
6. holm_recomputed        Holm adjustment across the four empirical P and the
                          four Fisher P per scenario matches stored values.
7. empirical_p_formula    empirical_p == (1 + exceedances) / (draws + 1).
8. flip_matrix_contract   the frozen donor sign-flip matrix has the registered
                          shape, all entries +/-1, and per-donor means ~ 0.
9. boundary_consistency   detectability boundary rows are internally consistent
                          (minimum >= observed; primary alt_up_pe_negative
                          minimum 3 over null mean 0.14856 > 20x).

A full replay of the 100,000-draw empirical null (re-running limma per flip)
was independently performed by the external module-3 audit and is NOT repeated
here; this script verifies everything downstream of the frozen flip matrix.

Read-only with respect to all frozen directories; the report is written to
work/03_bidirectional_convergence/independent_numeric_verification_v6_20260910/.

Usage:
    python3 scripts/03_bidirectional_convergence/verify_v6_numeric_independent.py
        [--project-root DIR] [--report-dir DIR]
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

ANALYSIS_ID = "m3_exploratory_shared_programs_v6_20260910"
RELEASE_REL = "data/processed/03_bidirectional_convergence/releases/" + ANALYSIS_ID
RESULTS_REL = "results/03_bidirectional_convergence/analyses/" + ANALYSIS_ID
TOL = 1e-9


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def close(a, b, tol=TOL):
    x, y = float(a), float(b)
    return abs(x - y) <= tol * max(1.0, abs(x), abs(y))


def hypergeom_sf(k, n_total, n_success, n_draw):
    """P(X >= k) for X ~ Hypergeometric(n_total, n_success, n_draw)."""
    if k <= 0:
        return 1.0
    lo = max(0, n_draw + n_success - n_total)
    hi = min(n_draw, n_success)
    if k > hi:
        return 0.0
    denom = math.comb(n_total, n_draw)
    return sum(math.comb(n_success, i) * math.comb(n_total - n_success, n_draw - i)
               for i in range(max(k, lo), hi + 1)) / denom


def holm(pvals):
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    m = len(pvals)
    out = [0.0] * m
    prev = 0.0
    for rank, i in enumerate(order):
        adj = (m - rank) * pvals[i]
        prev = max(prev, adj)
        out[i] = min(1.0, prev)
    return out


class Verifier:
    def __init__(self):
        self.checks = {}
        self.failures = []

    def record(self, name, ok, detail=""):
        self.checks[name] = bool(ok)
        if not ok:
            self.failures.append(f"{name}: {detail}" if detail else name)
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail and not ok else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    ap.add_argument("--project-root", default=str(root))
    ap.add_argument("--report-dir", default=None)
    args = ap.parse_args()
    root = Path(args.project_root).resolve()
    release = root / RELEASE_REL
    results = root / RESULTS_REL
    ver = Verifier()

    registry = json.loads((release / "module3_v5_analysis_registry.json").read_text())

    # 1. input hashes
    hash_problems = []
    for group in ("input_sources", "locked_outputs"):
        for name, meta in registry[group].items():
            p = Path(meta["path"])
            if not p.is_file():
                hash_problems.append(f"{group}/{name}: missing {p}")
            elif sha256_of(p) != meta["sha256"]:
                hash_problems.append(f"{group}/{name}: sha256 mismatch")
    ver.record("input_hashes", not hash_problems, "; ".join(hash_problems[:5]))

    # 2. universe contract
    universe = read_tsv(release / "gene_universe_and_mapping.tsv")
    u_problems = []
    if len(universe) != int(registry["universe_n"]) or len(universe) != 8101:
        u_problems.append(f"universe rows {len(universe)} != 8101")
    sym = [r["hgnc_symbol"] for r in universe]
    hid = [r["hgnc_id"] for r in universe]
    hg = [r["hgnc_ensembl_gene_id"] for r in universe]
    jti = [r["jti_ensembl_gene_id"] for r in universe]
    if len(set(sym)) != len(sym):
        u_problems.append("duplicate hgnc_symbol in universe")
    if len(set(hid)) != len(hid):
        u_problems.append("duplicate hgnc_id in universe")
    if len(set(hg)) != len(hg) or len(set(jti)) != len(jti) or set(hg) != set(jti):
        u_problems.append("Ensembl mapping not one-to-one across PE(JTI)/GSE103927 sides")
    ver.record("universe_contract", not u_problems, "; ".join(u_problems))

    # 3. strict sets rebuilt from raw sources
    alt_stats = read_tsv(Path(registry["input_sources"]["formal_alt_statistics"]["path"]))
    alt_by_symbol = {}
    dup = 0
    for r in alt_stats:
        s = r["hgnc_symbol"]
        if s in alt_by_symbol:
            dup += 1
        alt_by_symbol[s] = r
    pe_stats = read_tsv(Path(registry["input_sources"]["pe_whole_blood"]["path"]))
    pe_by_gene = {r["gene_id"]: r for r in pe_stats}

    u_syms = set(sym)
    alt_up = {r["hgnc_symbol"] for r in universe
              if r["hgnc_symbol"] in alt_by_symbol
              and float(alt_by_symbol[r["hgnc_symbol"]]["module2_fdr"]) < 0.05
              and float(alt_by_symbol[r["hgnc_symbol"]]["module2_moderated_logFC"]) >= 0.5}
    alt_down = {r["hgnc_symbol"] for r in universe
                if r["hgnc_symbol"] in alt_by_symbol
                and float(alt_by_symbol[r["hgnc_symbol"]]["module2_fdr"]) < 0.05
                and float(alt_by_symbol[r["hgnc_symbol"]]["module2_moderated_logFC"]) <= -0.5}
    pe_pos = {r["hgnc_symbol"] for r in universe
              if r["jti_ensembl_gene_id"] in pe_by_gene
              and float(pe_by_gene[r["jti_ensembl_gene_id"]]["q_gene_tissue"]) < 0.05
              and float(pe_by_gene[r["jti_ensembl_gene_id"]]["zscore"]) > 0}
    pe_neg = {r["hgnc_symbol"] for r in universe
              if r["jti_ensembl_gene_id"] in pe_by_gene
              and float(pe_by_gene[r["jti_ensembl_gene_id"]]["q_gene_tissue"]) < 0.05
              and float(pe_by_gene[r["jti_ensembl_gene_id"]]["zscore"]) < 0}

    overlap = read_tsv(results / "directional_overlap.tsv")
    primary = {r["quadrant"]: r for r in overlap if r["scenario"] == "PRIMARY_STRICT"}
    set_problems = []
    if dup:
        set_problems.append(f"{dup} duplicate symbols in ALT statistics (join ambiguous)")
    rebuilt = {
        "alt_up_pe_positive": (alt_up, pe_pos, len(alt_up), 9),
        "alt_up_pe_negative": (alt_up, pe_neg, len(alt_up), 91),
        "alt_down_pe_positive": (alt_down, pe_pos, len(alt_down), 35),
        "alt_down_pe_negative": (alt_down, pe_neg, len(alt_down), 12),
    }
    # fix expected ordering: quadrant -> (alt_n, pe_n)
    expected_n = {"alt_up_pe_positive": (91, 9), "alt_up_pe_negative": (91, 12),
                  "alt_down_pe_positive": (35, 9), "alt_down_pe_negative": (35, 12)}
    for q, (a, p) in expected_n.items():
        got = rebuilt[q]
        if len(got[0]) != a or len(got[1]) != p:
            set_problems.append(f"{q}: rebuilt sets {len(got[0])}/{len(got[1])} != expected {a}/{p}")
        if int(primary[q]["altitude_set_n"]) != len(got[0]) or int(primary[q]["pe_set_n"]) != len(got[1]):
            set_problems.append(f"{q}: stored set sizes disagree with rebuilt")
    ver.record("strict_sets_rebuilt", not set_problems, "; ".join(set_problems[:5]))

    # 4. observed intersections rebuilt
    inter_problems = []
    obs = {"alt_up_pe_positive": len(alt_up & pe_pos), "alt_up_pe_negative": len(alt_up & pe_neg),
           "alt_down_pe_positive": len(alt_down & pe_pos), "alt_down_pe_negative": len(alt_down & pe_neg)}
    for q, k in obs.items():
        if int(primary[q]["observed_overlap_n"]) != k:
            inter_problems.append(f"{q}: stored {primary[q]['observed_overlap_n']} != rebuilt {k}")
    wars1 = (alt_up & pe_neg)
    if wars1 != {"WARS1"}:
        inter_problems.append(f"alt_up_pe_negative intersection != {{WARS1}}: {sorted(wars1)}")
    genes = read_tsv(results / "directional_overlap_genes.tsv")
    prim_genes = {r["quadrant"] for r in genes if r["scenario"] == "PRIMARY_STRICT"}
    if prim_genes != {"alt_up_pe_negative"} or any(r["hgnc_symbol"] != "WARS1" for r in genes if r["scenario"] == "PRIMARY_STRICT"):
        inter_problems.append("PRIMARY_STRICT gene table does not record exactly WARS1 in alt_up_pe_negative")
    ver.record("observed_intersections", not inter_problems, "; ".join(inter_problems[:5]))

    # 5. Fisher one-sided recomputed for all 12 rows
    f_problems = []
    for r in overlap:
        n = int(r["universe_n"]); a = int(r["altitude_set_n"]); p = int(r["pe_set_n"])
        k = int(r["observed_overlap_n"])
        pval = hypergeom_sf(k, n, p, a)
        if not close(pval, r["fisher_p_one_sided"]):
            f_problems.append(f"{r['scenario']}/{r['quadrant']}: fisher {r['fisher_p_one_sided']} != {pval}")
    ver.record("fisher_recomputed", not f_problems, "; ".join(f_problems[:5]))

    # 6+7. Holm recompute + empirical P formula, per scenario family
    h_problems = []
    scenarios = sorted({r["scenario"] for r in overlap})
    for sc in scenarios:
        rows = [r for r in overlap if r["scenario"] == sc]
        emp = [float(r["empirical_p_one_sided"]) for r in rows]
        fis = [float(r["fisher_p_one_sided"]) for r in rows]
        h_emp, h_fis = holm(emp), holm(fis)
        for r, he, hf in zip(rows, h_emp, h_fis):
            if not close(he, r["empirical_holm_p"]):
                h_problems.append(f"{sc}/{r['quadrant']}: empirical Holm {r['empirical_holm_p']} != {he}")
            if not close(hf, r["fisher_holm_p"]):
                h_problems.append(f"{sc}/{r['quadrant']}: Fisher Holm {r['fisher_holm_p']} != {hf}")
            draws = int(r["empirical_draws"]); exc = int(r["empirical_exceedances"])
            if not close((1 + exc) / (draws + 1), r["empirical_p_one_sided"]):
                h_problems.append(f"{sc}/{r['quadrant']}: empirical P != (1+exc)/(draws+1)")
            if draws != int(registry["sign_flip_draws"]):
                h_problems.append(f"{sc}/{r['quadrant']}: draws {draws} != registry")
    ver.record("holm_and_empirical_p_formula", not h_problems, "; ".join(h_problems[:5]))

    # 8. flip matrix contract
    import numpy as np
    m_problems = []
    npz_path = Path(registry["locked_outputs"]["donor_sign_flip_matrix"]["path"])
    with np.load(npz_path) as z:
        keys = list(z.keys())
        mat = z[keys[0]]
    donors, draws = int(registry["donor_n"]), int(registry["sign_flip_draws"])
    if mat.shape not in ((draws, donors), (donors, draws)):
        m_problems.append(f"flip matrix shape {mat.shape} not ({draws},{donors}) or transpose")
    if not np.all((mat == 1) | (mat == -1)):
        m_problems.append("flip matrix has entries other than +/-1")
    if np.abs(mat.mean(axis=0)).max() > 0.02 if mat.shape[0] == draws else np.abs(mat.mean(axis=1)).max() > 0.02:
        m_problems.append("some donor's flip mean deviates from 0 by > 0.02")
    ver.record("flip_matrix_contract", not m_problems, "; ".join(m_problems[:5]))

    # 9. boundary consistency
    boundary = read_tsv(results / "detectability_boundary.tsv")
    b_problems = []
    for r in boundary:
        if int(r["minimum_overlap_for_bonferroni_upper_bound_lt_0_05"]) < int(r["observed_overlap_n"]):
            b_problems.append(f"{r['scenario']}/{r['quadrant']}: boundary below observed overlap")
    prim_neg = next(r for r in boundary if r["scenario"] == "PRIMARY_STRICT" and r["quadrant"] == "alt_up_pe_negative")
    ratio = int(prim_neg["minimum_overlap_for_bonferroni_upper_bound_lt_0_05"]) / float(prim_neg["null_mean_overlap"])
    if not (close(float(prim_neg["null_mean_overlap"]), 0.14856, 1e-4)
            and int(prim_neg["minimum_overlap_for_bonferroni_upper_bound_lt_0_05"]) == 3 and ratio > 20):
        b_problems.append(f"alt_up_pe_negative boundary inconsistent (min={prim_neg['minimum_overlap_for_bonferroni_upper_bound_lt_0_05']}, ratio={ratio:.2f})")
    ver.record("boundary_consistency", not b_problems, "; ".join(b_problems[:5]))

    report = {
        "analysis_id": ANALYSIS_ID,
        "purpose": "audit B-12: independent numeric verification of frozen module-3 v6 results",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "checks": ver.checks,
        "checks_passed": sum(ver.checks.values()),
        "checks_total": len(ver.checks),
        "rebuilt_primary_sets": {q: {"altitude": n[0], "pe": n[1]} for q, n in expected_n.items()},
        "rebuilt_primary_intersections": obs,
        "empirical_null_replay": "not repeated here; independently verified by the external module-3 audit (2026-09-10)",
        "failures": ver.failures,
        "status": "PASS" if not ver.failures else "FAIL",
    }
    report_dir = Path(args.report_dir).resolve() if args.report_dir else (
        root / "work/03_bidirectional_convergence/independent_numeric_verification_v6_20260910")
    report_dir.mkdir(parents=True, exist_ok=True)
    out = report_dir / "report.json"
    out.write_text(json.dumps(report, indent=1) + "\n")
    print(f"{report['status']}: {report['checks_passed']}/{report['checks_total']} checks; report: {out}")
    return 0 if not ver.failures else 1


if __name__ == "__main__":
    sys.exit(main())
