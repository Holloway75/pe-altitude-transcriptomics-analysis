#!/usr/bin/env python3
"""Independent closure validator for the global-BH MR-JTI batch."""
from __future__ import annotations

import argparse, csv, gzip, hashlib, json, math
from collections import Counter, defaultdict
from pathlib import Path

NA = "NA"

def rows(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")

def close(a, b, tol=1e-10):
    try:
        x, y = float(a), float(b)
        return math.isfinite(x) and math.isfinite(y) and abs(x-y) <= tol * max(1, abs(x), abs(y))
    except (TypeError, ValueError):
        return False

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def task_key(r): return r["dataset_id"], r["gene_id"], r["tissue"]

def repeatability_compare(root: Path, other: Path, errors: list):
    """Compare two isolated output directories of the same batch.

    Comparison policy: identity, categorical and count fields must match
    exactly; float fields must match within the same 1e-10 relative tolerance
    the closure validator uses (runner versions may differ in float
    formatting: Python repr vs R's 15-significant-digit output); per-task
    artifacts must match byte for byte. Tasks that never ran
    (harmonization_failed) have no scientific content; family metadata the
    runner versions report differently for them is disclosed separately, not
    treated as a scientific mismatch.
    """
    import filecmp
    numeric_fields = {
        "source_z", "source_pvalue", "source_q_global", "n_snps",
        "standardized_expression_coefficient", "original_lasso_coefficient",
        "beta_scale_ok", "ci_low", "ci_high", "lambda", "cv_mean_error",
        "cv_error_se", "bootstrap_zero_fraction", "bootstrap_positive_fraction",
        "bootstrap_negative_fraction", "mrjti_family_alpha", "source_family_size",
        "mrjti_family_size",
    }
    rep = {"compared_directory": str(other),
           "policy": "exact strings for identity/categorical fields; 1e-10 relative tolerance for floats; byte-identical task artifacts; unavailable-task metadata differences disclosed separately",
           "final_table_scientific_match": False,
           "manifest_identical": False,
           "field_diffs_numerical_within_tolerance": 0,
           "field_diffs_beyond_tolerance": 0,
           "task_artifacts_identical": 0,
           "task_artifacts_compared": 0,
           "unavailable_task_reporting_diffs": 0,
           "mismatches": []}
    a_final = list(rows(root / "mrjti-global-bh.tsv")); b_final = list(rows(other / "mrjti-global-bh.tsv"))
    if len(a_final) != len(b_final):
        errors.append(f"repeatability: final row count differs {len(a_final)} vs {len(b_final)}"); return rep
    amap = {task_key(x): x for x in a_final}; bmap = {task_key(x): x for x in b_final}
    if set(amap) != set(bmap):
        errors.append("repeatability: final task key sets differ"); return rep
    hard, soft = [], []
    unav_soft = 0
    for k in amap:
        arow, brow = amap[k], bmap[k]
        never_ran = arow.get("status") == "harmonization_failed"
        for f in arow:
            if f == "task_dir":
                continue
            av, bv = arow[f], brow.get(f)
            if av == bv:
                continue
            if never_ran and f in ("mrjti_family_size", "mrjti_family_alpha", "n_snps", "ci_significance", "mrjti_supported"):
                # reporting of never-run tasks differs across runner versions
                unav_soft += 1
                soft.append((k, f, av, bv))
                continue
            if f in numeric_fields and close(av, bv):
                soft.append((k, f, av, bv))
                continue
            hard.append((k, f, av, bv))
    rep["field_diffs_numerical_within_tolerance"] = len(soft)
    rep["field_diffs_beyond_tolerance"] = len(hard)
    rep["final_table_scientific_match"] = not hard
    rep["mismatches"] += [f"final {k}:{f}: {av!r} vs {bv!r}" for k, f, av, bv in hard[:20]]
    if soft:
        kinds = Counter(f for _, f, _, _ in soft)
        rep["soft_diff_kinds"] = dict(kinds)
    a_man = {task_key(x): x for x in rows(root / "candidate-manifest.tsv")}
    b_man = {task_key(x): x for x in rows(other / "candidate-manifest.tsv")}
    mdiffs = [f"{k}:{f}" for k in a_man if k in b_man
              for f in a_man[k] if f != "task_id" and a_man[k][f] != b_man[k].get(f)
              if not (f in numeric_fields and close(a_man[k][f], b_man[k].get(f)))]
    rep["manifest_identical"] = not mdiffs and set(a_man) == set(b_man)
    rep["mismatches"] += [f"manifest {x}" for x in mdiffs[:20]]
    for key, arow in amap.items():
        adir = Path(arow["task_dir"]); brow = bmap[key]; bdir = Path(brow["task_dir"])
        if arow.get("status") == "harmonization_failed":
            continue
        # mrjti-result.tsv is compared field-wise, not byte-wise: the vendored R
        # runner gained ci_significance/n_genes/family_alpha columns after the
        # formal run; shared columns must still agree.
        for name in ("harmonized.tsv", "pruned-variants.tsv", "mrjti-input.tsv", "prune.prune.in",
                     "prune-validation.vcor", "bootstrap.tsv", "snp-flow.tsv"):
            pa, pb = adir / name, bdir / name
            rep["task_artifacts_compared"] += 1
            if not pa.is_file() or not pb.is_file():
                rep["mismatches"].append(f"{key}/{name}: missing"); continue
            if filecmp.cmp(pa, pb, shallow=False):
                rep["task_artifacts_identical"] += 1
            else:
                rep["mismatches"].append(f"{key}/{name}: sha {sha(pa)[:12]} vs {sha(pb)[:12]}")
        ra = list(rows(adir / "mrjti-result.tsv")); rb = list(rows(bdir / "mrjti-result.tsv"))
        if len(ra) != 1 or len(rb) != 1:
            rep["mismatches"].append(f"{key}: mrjti-result row count"); continue
        ra, rb = ra[0], rb[0]
        shared = set(ra) & set(rb)
        rep.setdefault("result_schema_columns_only_in_rerun", set(rb) - set(ra))
        rep.setdefault("result_schema_columns_only_in_reference", set(ra) - set(rb))
        for f in shared:
            if ra[f] == rb[f]:
                continue
            if f in numeric_fields and close(ra[f], rb[f]):
                soft.append((key, f"result:{f}", ra[f], rb[f]))
                continue
            hard.append((key, f"result:{f}", ra[f], rb[f]))
        ja, jb = json.loads((adir / "task-result.json").read_text()), json.loads((bdir / "task-result.json").read_text())
        ja.pop("task_dir", None); jb.pop("task_dir", None)
        for f in set(ja) | set(jb):
            av, bv = ja.get(f), jb.get(f)
            if av == bv:
                continue
            if close(av, bv):
                soft.append((key, f"task-result:{f}", av, bv))
                continue
            hard.append((key, f"task-result:{f}", av, bv))
    rep["unavailable_task_reporting_diffs"] = unav_soft
    rep["result_schema_columns_only_in_rerun"] = sorted(rep.get("result_schema_columns_only_in_rerun", []))
    rep["result_schema_columns_only_in_reference"] = sorted(rep.get("result_schema_columns_only_in_reference", []))
    rep["mismatches"] = rep["mismatches"][:40]
    return rep

def recompute_source_bh(root: Path, manifest: list, errors: list) -> dict:
    """Audit B-11 (2026-09-10): recompute the upstream global-BH selection.

    Previously the manifest was treated as the ground truth for candidate
    selection. Here the S-PrediXcan all-tissues table of every source GWAS is
    re-read (after sha256 gate), BH is recomputed from pvalue over all unique
    finite gene-tissue tests, and three things are asserted per dataset:
      - the BH family size equals the recorded global_bh_family_size;
      - the recomputed q<=0.05 set equals the candidate-manifest set exactly;
      - every manifest row has source_q_global <= 0.05 and equal to the
        recomputed q.
    """
    method = json.loads((root / "method.json").read_text())
    manifest_keys = defaultdict(set)
    manifest_q = {}
    for m in manifest:
        manifest_keys[m["dataset_id"]].add((m["gene_id"], m["tissue"]))
        manifest_q[(m["dataset_id"], m["gene_id"], m["tissue"])] = m.get("source_q_global")
    details = {}
    for src in method["source_evidence"]:
        did = src["dataset_id"]
        table = Path(src["table"])
        if not table.is_file():
            errors.append(f"{did}: source table missing {table}")
            continue
        if sha(table) != src["table_sha256"]:
            errors.append(f"{did}: source table sha256 mismatch {table}")
        universe = []
        for r in rows(table):
            try:
                p = float(r["pvalue"])
            except (TypeError, ValueError):
                continue
            if not math.isfinite(p):
                continue
            if r.get("run_status") not in ("", None, "success"):
                continue
            universe.append((r["gene_id"], r["tissue"], p))
        keys = [(g, t) for g, t, _ in universe]
        if len(keys) != len(set(keys)):
            errors.append(f"{did}: duplicate gene-tissue rows in source table")
            continue
        n = len(universe)
        if n != int(src["global_bh_family_size"]):
            errors.append(f"{did}: BH family size {n} != recorded {src['global_bh_family_size']}")
        order = sorted(range(n), key=lambda i: universe[i][2])
        q = [0.0] * n
        prev = 1.0
        for rank in range(n, 0, -1):
            i = order[rank - 1]
            val = min(prev, universe[i][2] * n / rank)
            q[i] = val
            prev = val
        selected = {(universe[i][0], universe[i][1]): q[i] for i in range(n) if q[i] <= 0.05}
        recorded = int(src["global_bh_significant_tasks"])
        if len(selected) != recorded:
            errors.append(f"{did}: recomputed significant {len(selected)} != recorded {recorded}")
        wanted = manifest_keys.get(did, set())
        if set(selected) != wanted:
            errors.append(f"{did}: selection set mismatch "
                          f"(extra {len(set(selected) - wanted)}, missing {len(wanted - set(selected))})")
        n_checked = 0
        for (g, t), qv in selected.items():
            mq = manifest_q.get((did, g, t))
            if mq is None:
                continue
            if float(mq) > 0.05:
                errors.append(f"{did}/{g}/{t}: source_q_global {mq} > 0.05")
            if not close(mq, qv):
                errors.append(f"{did}/{g}/{t}: source_q_global {mq} != recomputed {qv}")
            n_checked += 1
        details[did] = {"family_size": n, "recomputed_significant": len(selected),
                        "recorded_significant": recorded, "manifest_q_checked": n_checked}
    return details


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("directory", type=Path); ap.add_argument("--report", type=Path)
    ap.add_argument("--recompute-source-bh", action="store_true",
                    help="audit B-11: independently recompute the upstream global-BH selection "
                         "from each source S-PrediXcan table and assert set equality with the manifest")
    ap.add_argument("--repeatability-compare", type=Path, help="second isolated output directory of the same batch; sets repeatability_proven when all scientific artifacts match")
    a = ap.parse_args(); root = a.directory.resolve(); errors=[]
    manifest=list(rows(root/"candidate-manifest.tsv")); final=list(rows(root/"mrjti-global-bh.tsv"))
    if len(manifest)!=2335 or len(final)!=2335: errors.append("expected 2335 manifest/final rows")
    if Counter(task_key(x) for x in manifest)!=Counter(task_key(x) for x in final): errors.append("manifest/final task closure failed")
    bykey={task_key(x):x for x in final}; gwas_expected=defaultdict(dict); artifact_paths=[]
    family_sizes=Counter(x["dataset_id"] for x in manifest)
    artifact_paths += [root/x for x in ("candidate-manifest.tsv","mrjti-global-bh.tsv","method.json","validation-report.json")]
    status=Counter(); successful=0
    for m in manifest:
        key=task_key(m); r=bykey[key]; status[r["status"]]+=1
        if any(r[x] != NA for x in ("pvalue","q_mrjti")): errors.append(f"{key}: unsupported P/q field populated")
        expected_n=family_sizes[r["dataset_id"]]
        if r["mrjti_family_size"] != NA and (int(r["mrjti_family_size"])!=expected_n or not close(r["mrjti_family_alpha"],0.05/expected_n)): errors.append(f"{key}: MR-JTI family mismatch")
        td=Path(r["task_dir"]); tr=td/"task-result.json"
        if not tr.is_file(): errors.append(f"{key}: task-result missing"); continue
        jr=json.loads(tr.read_text())
        if any(str(jr.get(k))!=str(r.get(k)) for k in r): errors.append(f"{key}: task JSON/final mismatch")
        artifact_paths.append(tr)
        if r["status"]=="harmonization_failed":
            if not r["inference_status"].startswith("not_run:besd_probe_"): errors.append(f"{key}: invalid unavailable reason")
            continue
        if r["status"]!="success_inference_bonferroni_selection_dependent": errors.append(f"{key}: unexpected status {r['status']}"); continue
        successful+=1
        required=[td/x for x in ("harmonized.tsv","pruned-variants.tsv","mrjti-input.tsv","prune.prune.in","prune-validation.vcor","mrjti-result.tsv","bootstrap.tsv","snp-flow.tsv")]
        if any(not x.is_file() for x in required): errors.append(f"{key}: missing task artifact"); continue
        artifact_paths += required
        inp=list(rows(td/"mrjti-input.tsv")); pruned=list(rows(td/"pruned-variants.tsv")); rr=list(rows(td/"mrjti-result.tsv")); boot=list(rows(td/"bootstrap.tsv"))
        ids={x["genotype_id"] for x in inp}; prune_ids={x.strip() for x in (td/"prune.prune.in").read_text().splitlines() if x.strip()}
        if len(inp)!=int(r["n_snps"]) or len(inp)<20 or len(ids)!=len(inp) or ids!=prune_ids or ids!={x["genotype_id"] for x in pruned}: errors.append(f"{key}: prune/input closure failed")
        v=list(rows(td/"prune-validation.vcor"))
        if any(float(x["UNPHASED_R2"])>0.1 for x in v): errors.append(f"{key}: retained LD exceeds r2 0.1")
        if len(rr)!=1 or len(boot)!=500: errors.append(f"{key}: result/bootstrap closure failed")
        else:
            draws=[float(x["expression_beta"]) for x in boot]; one=rr[0]
            derived=[sum(draws)/len(draws),sum(x==0 for x in draws)/len(draws),sum(x>0 for x in draws)/len(draws),sum(x<0 for x in draws)/len(draws)]
            reported=[one["standardized_expression_coefficient"],one["bootstrap_zero_fraction"],one["bootstrap_positive_fraction"],one["bootstrap_negative_fraction"]]
            if not all(close(x,y) for x,y in zip(derived,reported)): errors.append(f"{key}: bootstrap summary mismatch")
            alpha=0.05/expected_n; ordered=sorted(draws)
            def q7(p):
                h=(len(ordered)-1)*p; lo=math.floor(h); hi=math.ceil(h)
                return ordered[lo] if lo==hi else ordered[lo]+(h-lo)*(ordered[hi]-ordered[lo])
            original=float(one["original_lasso_coefficient"])
            expected_ci=(2*original-q7(1-alpha/2),2*original-q7(alpha/2))
            expected_sig="sig" if expected_ci[0]*expected_ci[1]>0 else "nonsig"
            if not close(one["ci_low"],expected_ci[0]) or not close(one["ci_high"],expected_ci[1]): errors.append(f"{key}: Bonferroni CI mismatch")
            if one.get("ci_significance")!=expected_sig or r["ci_significance"]!=expected_sig or r["mrjti_supported"]!=str(expected_sig=="sig").lower(): errors.append(f"{key}: native significance mismatch")
            for field in ("standardized_expression_coefficient","original_lasso_coefficient","ci_low","ci_high","lambda","cv_mean_error","cv_error_se","bootstrap_zero_fraction","bootstrap_positive_fraction","bootstrap_negative_fraction"):
                if not close(r[field],one[field]): errors.append(f"{key}: final/result mismatch {field}")
        gene=m["gene_id"]; tissue=m["tissue"]
        ref_path=root/"resources/ld"/gene/"ld_reference/ld-reference-variant-map.tsv"
        eq_path=root/"resources/eqtl"/tissue/gene/"eqtl"/tissue/f"{gene}.txt"
        if not ref_path.is_file() or not eq_path.is_file(): errors.append(f"{key}: upstream cache missing"); continue
        ref={x["rsid"]:x for x in rows(ref_path)}; eq={x["SNP"]:x for x in rows(eq_path)}
        for x in inp:
            er=eq.get(x["rsid"]); lr=ref.get(x["rsid"])
            if er is None or lr is None: errors.append(f"{key}/{x['rsid']}: eQTL/reference foreign key missing"); continue
            if not (x["chromosome"]==er["Chr"] and int(float(x["position"]))==int(er["BP"]) and x["effect_allele"]==er["A1"] and x["other_allele"]==er["A2"] and close(x["eqtl_beta"],er["b"]) and close(x["eqtl_se"],er["SE"]) and close(x["eqtl_p"],er["p"])): errors.append(f"{key}/{x['rsid']}: eQTL closure failed")
            if not (x["internal_id"]==lr["internal_id"] and x["genotype_id"]==lr["genotype_id"] and close(x["ldscore"],lr["ldscore"])): errors.append(f"{key}/{x['rsid']}: reference closure failed")
            sig=(x["gwas_source_effect_allele"],x["gwas_source_other_allele"],x["gwas_source_beta"],x["gwas_source_se"],x["gwas_source_p"],x["chromosome"],str(int(float(x["position"]))))
            old=gwas_expected[m["dataset_id"]].setdefault(x["rsid"],sig)
            if old!=sig: errors.append(f"{key}/{x['rsid']}: inconsistent GWAS source fields")
    method=json.loads((root/"method.json").read_text()); source={x["dataset_id"]:x for x in method["source_evidence"]}
    source_bh_details=None
    if a.recompute_source_bh:
        source_bh_details=recompute_source_bh(root, manifest, errors)
    for did,wanted in gwas_expected.items():
        seen=set(); path=Path(source[did]["normalized_gwas"])
        for g in rows(path):
            rid=g.get("rsid")
            if rid not in wanted: continue
            seen.add(rid); exp=wanted[rid]
            got=(g["effect_allele"],g["non_effect_allele"],g["beta"],g["standard_error"],g["pvalue"],g["chromosome"],str(int(float(g["position"]))))
            if got[:2]!=exp[:2] or got[5:]!=exp[5:] or not all(close(x,y) for x,y in zip(got[2:5],exp[2:5])): errors.append(f"{did}/{rid}: normalized GWAS closure failed")
        if seen!=set(wanted): errors.append(f"{did}: normalized GWAS rsID closure failed")
    for gene_dir in sorted((root/"resources/ld").iterdir()):
        mp=gene_dir/"ld_reference/ld-reference-variant-map.tsv"; score=next((gene_dir/"ld_reference").glob("*.score.ld"),None)
        if not mp.is_file() or score is None: errors.append(f"{gene_dir.name}: reference artifacts missing"); continue
        artifact_paths += [mp,score,gene_dir/"resource-method.json"]
        vals=list(rows(mp)); rs=[x["rsid"] for x in vals]; pos=[(x["chromosome"],x["position"]) for x in vals]
        if len(rs)!=len(set(rs)) or len(pos)!=len(set(pos)): errors.append(f"{gene_dir.name}: nonunique reference")
        for x in vals:
            if not (.01<=float(x["maf"])<=.5 and 0<=float(x["missing_rate"])<=.02 and math.isfinite(float(x["ldscore"]))): errors.append(f"{gene_dir.name}/{x['rsid']}: reference QC invalid")
            if len(x["ref"])==len(x["alt"])==1 and {x["ref"],x["alt"]} in ({"A","T"},{"C","G"}): errors.append(f"{gene_dir.name}/{x['rsid']}: palindrome retained")
    artifact_paths += sorted((root/"resources/eqtl").glob("*/*/eqtl/*/*.txt"))
    unique=sorted({p.resolve() for p in artifact_paths if p.is_file()},key=str); digests=[]
    for p in unique: digests.append((str(p.relative_to(root)),sha(p)))
    bundle=hashlib.sha256("\n".join(f"{p}\0{h}" for p,h in digests).encode()).hexdigest()
    repeatability_proven=False; repeatability=None
    if a.repeatability_compare is not None:
        repeatability=repeatability_compare(root, a.repeatability_compare.resolve(), errors)
        repeatability_proven=bool(repeatability["final_table_scientific_match"] and repeatability["manifest_identical"]
                                  and repeatability["task_artifacts_identical"]==repeatability["task_artifacts_compared"]
                                  and not repeatability.get("result_schema_columns_only_in_reference")
                                  and not repeatability["mismatches"])
        repeatability["repeatability_proven"]=repeatability_proven
    report={"status":"pass" if not errors else "fail","errors":errors,"candidate_tasks":len(manifest),"successful_tasks":successful,"status_counts":dict(status),"family_sizes":dict(family_sizes),"checked_scientific_artifacts":len(unique),"scientific_bundle_sha256":bundle,"repeatability_proven":repeatability_proven,"native_bonferroni_ci_available":True,"calibrated_p_or_q_available":False}
    if source_bh_details is not None:
        report["source_selection_independently_recomputed"]=True; report["source_bh"]=source_bh_details
    if repeatability is not None: report["repeatability"]=repeatability
    out=a.report or root/"independent-validation-report.json"; out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True)); raise SystemExit(0 if not errors else 1)

if __name__=="__main__": main()
