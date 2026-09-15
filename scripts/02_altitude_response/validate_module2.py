#!/usr/bin/env python3
"""Read-only, value-level validator for the Module 2 scientific contract.

⚠ LEGACY / DO NOT USE FOR CURRENT DELIVERABLES (2026-09-10, audit B-10):
This validator implements a retired payload schema (SCHEMA_VERSION 5.0) whose
dataset whitelist ({GSE196728, GSE75665, GSE103940, GSE46480}) contains the four
datasets now EXCLUDED from the study per RESEARCH_PLAN §4.2, and it requires an
HGNC snapshot (hgnc_complete_set_20260727.tsv) that has been deleted from disk
(audit B-9; HGNC does not redistribute historical snapshots, so byte-level
re-freeze of that era is impossible — same recovery boundary as the module 3 v4
audit). It therefore cannot validate the current GSE103927/GSE333506 v2 axes
and is retained only as a historical record of the original contract.

Current deliverables are validated by validate_module2_v2_bh.py (independent BH
recompute + freeze manifest closure; report written to work/).
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math, re, sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Iterable

EXPECTED_DATASETS = {
    "GSE196728": {"tissue": "whole blood", "role": "WHOLE_BLOOD_DE"},
    "GSE75665": {"tissue": "whole blood", "role": "WHOLE_BLOOD_DE"},
    "GSE103940": {"tissue": "whole blood", "role": "WHOLE_BLOOD_DE"},
    "GSE46480": {"tissue": "pbmc", "role": "PBMC_VALIDATION"},
}
CONFIRMATION_CANDIDATES = {"GSE75665", "GSE103940"}
AUTHORIZED_DATASETS = set(EXPECTED_DATASETS)
EXPECTED_FLOW = {"GSE196728": (144, 8), "GSE75665": (20, 10), "GSE46480": (196, 98), "GSE103940": (22, 11)}
GSE196728_CLOCKS = ("02", "06", "10", "14", "18", "22")
TRUE, FALSE, MISSING = {"true", "yes", "1"}, {"false", "no", "0"}, {"", "na", "n/a", "nan", "null", "."}
HEX64 = re.compile(r"^[0-9a-f]{64}$")
META_REQ = {"dataset_id","sample_id","expression_column","donor_id","include","exclusion_reason","tissue","assay","altitude_label","pair_id","source_url","planned_imbalance_id","acetazolamide_raw"}
DESIGN_REQ = {"dataset_id","model_id","estimand","formula","contrast_id","contrast_definition","scientific_comparison_id","reference_level","n_samples","n_donors","design_rank","design_columns","residual_df","correlation_engine","correlation_structure","block_vector_path","block_vector_sha256","donor_key_sha256","random_effect_formula","primary","interpretation_boundary","registry_version","design_matrix_path","design_matrix_sha256","sample_order_sha256","contrast_vector_path","contrast_vector_sha256"}
DE_REQ = {"dataset_id","model_id","contrast_id","estimand","test_id","feature_id","feature_namespace","feature_namespace_version","gene_id","gene_symbol","mapping_status","log2fc","se","statistic","statistic_type","df","p_value","fdr","direction","deg","base_expression","n_samples","n_donors","annotation_version","method_version","test_status","not_tested_reason","bh_family_id","test_universe_sha256"}
SOURCE_REQ = {"dataset_id","source_kind","url","retrieved_at_utc","local_path","bytes","sha256","content_role","archive_check","access_status","record_status_at_utc","citation","terms_url","license_status","use_restrictions","privacy_handling","local_verification_status","remote_immutable_id","verification_reason","audit_exception_status","audit_exception_id","audit_exception_approved_by"}
QC_REQ = {"dataset_id","stage","metric","value","status","threshold","details"}
ELIG_REQ = {"dataset_id","evidence_layer","eligibility_status","n_samples_included","n_donors_included","pairing_complete","donor_independence_verified","tissue_match","primary_contrast_estimable","design_full_rank","residual_df_adequate","source_data_available","overlap_status","decision_reason","registry_version","analysis_role","expected_samples","expected_donors","minimum_donors","maximum_attrition","exclusion_count","exclusion_reasons_complete","de_authorized","no_de_reason","dataset_manifest_sha256","dataset_payload_sha256"}
LAYER_REQ = {"dataset_id","analysis_role","analysis_state","contrast_id","estimand","output_path","output_sha256","eligibility_status","registry_version"}
COMPOSITION_REQ = {"dataset_id","total_model_id","conditional_model_id","contrast_id","test_id","total_bulk_log2fc","total_bulk_se","total_bulk_statistic","total_bulk_fdr","conditional_log2fc","conditional_se","conditional_statistic","conditional_fdr","delta_log2fc","attenuation_fraction","direction_relation","sensitivity_class","source_de_sha256"}
PBMC_REQ = {"dataset_id","model_id","contrast_id","estimand","test_id","target_id","target_source_sha256","hypothesis_status","statistic","p_value","fdr","multiplicity_family_id","registry_version"}
HYPOTHESIS_REQ = {"target_id","registry_version","frozen_at_utc","source_dataset_id","source_model_id","source_contrast_id","source_estimand","source_test_id","source_de_sha256","hypothesis_family_id","pbmc_model_id","pbmc_contrast_id","pbmc_estimand","pbmc_test_id","multiplicity_family_id"}
META_EFFECT_REQ = {"dataset_id","source_family_id","model_id","contrast_id","scientific_comparison_id","test_id","estimand","effect","se","source_de_sha256"}
META_RESULT_REQ = {"meta_family_id","scientific_comparison_id","estimand","test_id","k_studies","meta_log2fc","meta_se","meta_statistic","df","p_value","fdr","q_heterogeneity","i2","tau2","direction_concordant","estimand_comparability_status","modified_hksj_low","modified_hksj_high","fixed_effect_log2fc","fixed_effect_se","source_effects_sha256"}
SCHEMA_VERSION = "5.0"

def schema_sha256():
    payload={"metadata":sorted(META_REQ),"design":sorted(DESIGN_REQ),"de":sorted(DE_REQ),"source":sorted(SOURCE_REQ),"eligibility":sorted(ELIG_REQ),"layer":sorted(LAYER_REQ),"composition":sorted(COMPOSITION_REQ),"pbmc":sorted(PBMC_REQ),"hypothesis":sorted(HYPOTHESIS_REQ),"meta_effect":sorted(META_EFFECT_REQ),"meta_result":sorted(META_RESULT_REQ)}
    return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(",",":")).encode()).hexdigest()

@dataclass(frozen=True)
class Diagnostic:
    severity: str; code: str; message: str; path: str = ""; row: int | None = None

class Report:
    def __init__(self, mode: str, generation_context=None):
        self.mode, self.diagnostics, self.metrics = mode, [], {}
        self.generated_at_utc=datetime.now(timezone.utc).isoformat()
        self.generation_context=generation_context
    def add(self, severity, code, message, path="", row=None): self.diagnostics.append(Diagnostic(severity,code,message,str(path),row))
    def error(self, code, message, path="", row=None): self.add("ERROR",code,message,path,row)
    def warning(self, code, message, path="", row=None): self.add("WARNING",code,message,path,row)
    @property
    def error_count(self): return sum(d.severity == "ERROR" for d in self.diagnostics)
    @property
    def warning_count(self): return sum(d.severity == "WARNING" for d in self.diagnostics)
    def payload(self):
        payload={"schema_version":SCHEMA_VERSION,"schema_sha256":schema_sha256(),"validator_sha256":sha(Path(__file__)),"mode":self.mode,"generated_at_utc":self.generated_at_utc,"status":"PASS" if not self.error_count else "FAIL","schema_pass":not self.error_count,"scientific_gate_pass":not self.error_count,"publishable":self.mode != "release" or (not self.error_count and bool(self.metrics.get("success_marker_verified"))),"phase_b_complete":self.metrics.get("phase_b_complete",False),"payload_root_sha256":self.metrics.get("payload_root_sha256"),"error_count":self.error_count,"warning_count":self.warning_count,"metrics":self.metrics,"diagnostics":[asdict(d) for d in self.diagnostics]}
        if self.generation_context is not None:payload["generation_context"]=self.generation_context
        return payload

def missing(v): return v is None or str(v).strip().lower() in MISSING
def boolean(v, field, path, row, report):
    x=str(v).strip().lower()
    if x in TRUE:return True
    if x in FALSE:return False
    report.error("INVALID_BOOLEAN",f"{field} has invalid boolean {v!r}",path,row); return None
def number(v, field, path, row, report, allow_na=False):
    if allow_na and missing(v): return None
    try: x=float(v)
    except (TypeError,ValueError): report.error("INVALID_NUMBER",f"{field} is not numeric: {v!r}",path,row); return None
    if not math.isfinite(x): report.error("NONFINITE_NUMBER",f"{field} must be finite",path,row); return None
    return x
def integer(v, field, path, row, report):
    x=number(v,field,path,row,report)
    if x is None:return None
    if x != int(x): report.error("INVALID_INTEGER",f"{field} must be integral",path,row); return None
    return int(x)
def read_tsv(path, report, required=None, allow_empty=False):
    if not path.is_file(): report.error("MISSING_FILE","Required TSV is absent",path); return [],[]
    try:
        with path.open(encoding="utf-8",newline="") as h:
            rd=csv.DictReader(h,delimiter="\t"); cols=rd.fieldnames or []; rows=[dict(r) for r in rd]
    except Exception as e: report.error("UNREADABLE_TSV",str(e),path); return [],[]
    if len(cols)!=len(set(cols)): report.error("DUPLICATE_COLUMNS","Duplicate header",path)
    if required and not required.issubset(cols): report.error("MISSING_COLUMNS",", ".join(sorted(required-set(cols))),path)
    if not rows and not allow_empty: report.error("EMPTY_TABLE","No data rows",path)
    return rows,cols
def load_json(path, report):
    try: obj=json.loads(path.read_text(encoding="utf-8"))
    except Exception as e: report.error("INVALID_JSON",str(e),path); return None
    if not isinstance(obj,dict): report.error("INVALID_JSON_ROOT","Expected object",path); return None
    return obj
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
def valid_hash(v): return bool(HEX64.fullmatch(str(v))) and set(str(v)) != {"0"}
def valid_utc(v):
    try:
        text=str(v); dt=datetime.fromisoformat(text.replace("Z","+00:00")); return dt.tzinfo is not None and dt.utcoffset().total_seconds()==0
    except (ValueError,TypeError,AttributeError): return False
def safe_child(root, text):
    p=PurePosixPath(text)
    if p.is_absolute() or ".." in p.parts or not text or "\\" in text:return None
    out=(root/p).resolve()
    try: out.relative_to(root.resolve())
    except ValueError:return None
    return out
def bh(ps):
    n=len(ps); out=[1.0]*n; run=1.0
    for j,i in reversed(list(enumerate(sorted(range(n),key=ps.__getitem__),1))):
        run=min(run,ps[i]*n/j); out[i]=min(1.0,run)
    return out

def medication01(value):
    x=str(value).strip().lower()
    if x in {"false","0","no"}: return 0
    if x in {"true","1","yes"}: return 1
    raise ValueError(value)

def renv_dependencies(record):
    """Return mandatory package names from standard renv dependency fields."""
    dependencies=set()
    if not isinstance(record,dict):return dependencies
    for field in ("Depends","Imports","LinkingTo"):
        value=record.get(field,[])
        if isinstance(value,str):parts=value.split(",")
        elif isinstance(value,list):parts=value
        elif value in (None,{}):parts=[]
        else:parts=[value]
        for part in parts:
            name=re.split(r"\s|\(",str(part).strip(),maxsplit=1)[0].strip()
            if name and name!="R":dependencies.add(name)
    return dependencies

def validate_metadata(path, report, dataset_id=None):
    rows,cols=read_tsv(path,report,META_REQ)
    row_ids={r.get("dataset_id","").strip() for r in rows if r.get("dataset_id","").strip()}
    if dataset_id is None:
        if len(row_ids)==1: dataset_id=next(iter(row_ids))
        else: report.error("AMBIGUOUS_METADATA_DATASET_ID",f"Expected one dataset_id, got {sorted(row_ids)}",path)
    if dataset_id and dataset_id not in AUTHORIZED_DATASETS:report.error("UNAUTHORIZED_DATASET_ID",dataset_id,path)
    included=[]; seen_s=set(); seen_e=set(); donors=defaultdict(list)
    for rn,r in enumerate(rows,2):
        if dataset_id and r.get("dataset_id","").strip()!=dataset_id: report.error("DATASET_ID_MISMATCH",dataset_id,path,rn)
        inc=boolean(r.get("include",""),"include",path,rn,report)
        if inc is False and missing(r.get("exclusion_reason")): report.error("MISSING_EXCLUSION_REASON","Excluded row needs reason",path,rn)
        if inc is not True: continue
        included.append((rn,r))
        for f in ("sample_id","expression_column","donor_id","tissue","assay","altitude_label","pair_id","source_url"):
            if missing(r.get(f)): report.error("MISSING_INCLUDED_VALUE",f,path,rn)
        for f,s in (("sample_id",seen_s),("expression_column",seen_e)):
            v=r.get(f,"").strip()
            if v in s: report.error("DUPLICATE_INCLUDED_KEY",f"{f}={v}",path,rn)
            s.add(v)
        donors[r.get("donor_id","").strip()].append((rn,r))
    if not included: report.error("NO_INCLUDED_SAMPLES","No included samples",path)
    if dataset_id in AUTHORIZED_DATASETS:
        observed={" ".join(r.get("tissue","").lower().replace("_"," ").split()) for _,r in included}
        if dataset_id=="GSE46480": tissue_ok=all(x=="pbmc" or "peripheral blood mononuclear" in x for x in observed)
        else: tissue_ok=observed=={"whole blood"}
        if not tissue_ok:report.error("TISSUE_ROLE_MISMATCH",f"{dataset_id}: {sorted(observed)}",path)
    for donor,items in donors.items():
        pairs={x[1].get("pair_id","").strip() for x in items}
        if pairs != {donor}: report.error("PAIR_ID_MISMATCH",f"{donor}: {sorted(pairs)}",path)
    if dataset_id == "GSE75665":
        for donor,items in donors.items():
            labels=[x[1]["altitude_label"].strip() for x in items]
            if sorted(labels)!=["5300m","plain"]: report.error("GSE75665_PAIR_GRID",f"{donor}: {labels}",path)
    elif dataset_id == "GSE46480":
        medication_states=[]
        for donor,items in donors.items():
            labels=[x[1]["altitude_label"].strip() for x in items]
            if sorted(labels)!=["baseline","day3"]: report.error("GSE46480_PRIMARY_GRID",f"{donor}: {labels}",path)
            vals=[]
            for rn,r in items:
                raw=r.get("acetazolamide_raw","")
                if missing(raw):
                    medication_states.append("missing")
                    continue
                medication_states.append("known")
                try: vals.append(medication01(raw))
                except ValueError: report.error("INVALID_ACETAZOLAMIDE","Use explicit false/true, 0/1, no/yes, or audited all-NA unavailable status",path,rn)
            if vals and len(set(vals))!=1: report.error("MEDICATION_WITHIN_DONOR_MISMATCH",donor,path)
        if set(medication_states)=={"missing","known"}:
            report.error("PARTIAL_ACETAZOLAMIDE_MAPPING","Medication labels must be complete or audited unavailable for every included donor",path)
        elif medication_states and set(medication_states)=={"missing"}:
            report.warning("ACETAZOLAMIDE_LABELS_UNAVAILABLE","All donor-level labels are NA; medication sensitivities must be NOT_RUN_DONOR_LABELS_NOT_PUBLIC",path)
    elif dataset_id == "GSE103940":
        for donor,items in donors.items():
            labels=[x[1]["altitude_label"].strip() for x in items]
            if sorted(labels)!=["highaltitude","plain"]:report.error("GSE103940_PAIR_GRID",f"{donor}: {labels}",path)
    elif dataset_id == "GSE196728":
        all_clocks={r.get("time_of_day","").strip() for _,r in included if not missing(r.get("time_of_day"))}
        if all_clocks!=set(GSE196728_CLOCKS):report.error("GSE196728_CLOCK_SET",f"Expected {list(GSE196728_CLOCKS)}, got {sorted(all_clocks)}",path)
        for donor,items in donors.items():
            cells=[]
            for rn,r in items:
                clock=r.get("time_of_day","").strip(); period=r.get("altitude_label","").strip()
                if missing(clock): report.error("TIME_OF_DAY_NOT_FROZEN",donor,path,rn)
                cells.append((period,clock))
            if len(cells)!=len(set(cells)): report.error("DUPLICATE_REPEATED_CELL",donor,path)
            periods={p for p,_ in cells}
            if periods != {"sea_level","3800m","5100m"}: report.error("GSE196728_PERIOD_GRID",f"{donor}: {sorted(periods)}",path)
            expected={(p,c) for p in ("sea_level","3800m","5100m") for c in GSE196728_CLOCKS}
            missing_cells=expected-set(cells)
            if missing_cells and all(missing(r.get("planned_imbalance_id")) for _,r in items): report.error("UNDECLARED_IMBALANCE",f"{donor}: {sorted(missing_cells)}",path)
    return {"rows":rows,"included":included,"n_samples":len(included),"n_donors":len(donors),"donors":donors}

def matrix_rank(a,tol=1e-10):
    a=[list(map(float,row)) for row in a]; m=len(a); n=len(a[0]) if m else 0; rank=0
    for col in range(n):
        pivot=max(range(rank,m),key=lambda i:abs(a[i][col]),default=rank)
        if rank>=m or abs(a[pivot][col])<=tol: continue
        a[rank],a[pivot]=a[pivot],a[rank]; q=a[rank][col]
        for j in range(col,n):a[rank][j]/=q
        for i in range(m):
            if i!=rank:
                q=a[i][col]
                for j in range(col,n):a[i][j]-=q*a[rank][j]
        rank+=1
    return rank
def read_matrix(path,report):
    rows,cols=read_tsv(path,report)
    if not rows or len(cols)<2:return [],[],[]
    labels=[]; vals=[]
    for rn,r in enumerate(rows,2):
        labels.append(r[cols[0]])
        row=[]
        for c in cols[1:]:
            x=number(r[c],c,path,rn,report); row.append(0.0 if x is None else x)
        vals.append(row)
    return labels,cols[1:],vals
def validate_design(path,report,dataset_id,meta):
    rows,cols=read_tsv(path,report,DESIGN_REQ); keys={}
    sample_order=[r["sample_id"].strip() for _,r in meta["included"]]
    sample_hash=hashlib.sha256(("\n".join(sample_order)+"\n").encode()).hexdigest()
    donor_order=[r["donor_id"].strip() for _,r in meta["included"]]
    donor_hash=hashlib.sha256(("\n".join(donor_order)+"\n").encode()).hexdigest()
    for rn,r in enumerate(rows,2):
        if r.get("dataset_id")!=dataset_id:report.error("DESIGN_DATASET_ID_MISMATCH",f"Expected {dataset_id}",path,rn)
        key=(r.get("model_id","").strip(),r.get("contrast_id","").strip(),r.get("estimand","").strip())
        if key in keys:report.error("DUPLICATE_DESIGN_KEY",repr(key),path,rn)
        else:keys[key]=r
        for f in ("model_id","contrast_id","estimand","formula","scientific_comparison_id","registry_version"):
            if missing(r.get(f)):report.error("EMPTY_DESIGN_FIELD",f,path,rn)
        for f,actual in (("n_samples",meta["n_samples"]),("n_donors",meta["n_donors"])):
            x=integer(r.get(f,""),f,path,rn,report)
            if x is not None and x!=actual:report.error("COUNT_MISMATCH",f,path,rn)
        mp=safe_child(path.parent,r.get("design_matrix_path","")); cp=safe_child(path.parent,r.get("contrast_vector_path",""))
        if not mp or not cp: report.error("UNSAFE_EVIDENCE_PATH","matrix/contrast path",path,rn); continue
        for p,f in ((mp,"design_matrix_sha256"),(cp,"contrast_vector_sha256")):
            if not p.is_file() or sha(p)!=r.get(f):report.error("EVIDENCE_HASH_MISMATCH",f,path,rn)
        labels,names,x=read_matrix(mp,report)
        _,cnames,c=read_matrix(cp,report)
        rank=matrix_rank(x); n=len(x); p=len(names); rdf=n-rank
        if labels!=sample_order:report.error("DESIGN_SAMPLE_ORDER_MISMATCH","matrix rows",mp)
        if r.get("sample_order_sha256")!=sample_hash:report.error("SAMPLE_ORDER_HASH_MISMATCH","design",path,rn)
        for f,val in (("design_rank",rank),("design_columns",p),("residual_df",rdf)):
            z=integer(r.get(f,""),f,path,rn,report)
            if z is not None and z!=val:report.error("RECOMPUTED_DESIGN_MISMATCH",f,path,rn)
        if rank!=p:report.error("RANK_DEFICIENT_DESIGN","rank < columns",path,rn)
        if rdf<=0:report.error("NO_RESIDUAL_DF","residual df <= 0",path,rn)
        if dataset_id=="GSE75665" and ("subject:altitude" in r.get("formula","").replace(" ","") or p>=n):
            report.error("GSE75665_SATURATED_DESIGN","Use a positive-residual-df paired additive model, not donor-by-altitude cell means",path,rn)
        donors=donor_order;levels=sorted(set(donors)); donor_indicators=[[int(d==level) for level in levels[1:]] for d in donors]
        augmented=[row+ind for row,ind in zip(x,donor_indicators)]
        has_donor_fixed=matrix_rank(augmented)==rank
        engine=r.get("correlation_engine","")
        if r.get("donor_key_sha256")!=donor_hash:report.error("DONOR_KEY_HASH_MISMATCH","",path,rn)
        if engine=="fixed_effect":
            if not has_donor_fixed:report.error("DONOR_BLOCK_NOT_IN_DESIGN","Fixed-effect branch needs donor indicators",mp)
            if not missing(r.get("block_vector_path")) or not missing(r.get("random_effect_formula")):report.error("CONFLICTING_CORRELATION_EVIDENCE","Fixed effect cannot also declare block/random effect",path,rn)
        elif engine in {"duplicateCorrelation","dream"}:
            if has_donor_fixed:report.error("CONFLICTING_CORRELATION_EVIDENCE","Correlated/random branch must not also include full donor indicators",path,rn)
            bp=safe_child(path.parent,r.get("block_vector_path",""))
            if not bp or not bp.is_file() or sha(bp)!=r.get("block_vector_sha256"):report.error("BLOCK_VECTOR_HASH_MISMATCH","",path,rn)
            else:
                br,_=read_tsv(bp,report,{"sample_id","donor_id"})
                if [x.get("sample_id","") for x in br]!=sample_order or [x.get("donor_id","") for x in br]!=donor_order:report.error("BLOCK_VECTOR_MISMATCH","",bp)
            if engine=="duplicateCorrelation" and not missing(r.get("random_effect_formula")):report.error("CONFLICTING_CORRELATION_EVIDENCE","duplicateCorrelation does not use random-effect formula",path,rn)
            if engine=="dream" and r.get("random_effect_formula","").replace(" ","")!="(1|donor_id)":report.error("INVALID_RANDOM_EFFECT_FORMULA","",path,rn)
        else:report.error("INVALID_CORRELATION_ENGINE",engine,path,rn)
        if not c or len(c)!=1 or cnames!=names:report.error("INVALID_CONTRAST_VECTOR","Contrast columns must match design",cp)
        elif all(abs(v)<1e-12 for v in c[0]):report.error("ZERO_CONTRAST","Contrast is zero",cp)
    return keys

def validate_method(path,report,dataset_id):
    m=load_json(path,report)
    if not m:return {}
    required={"schema_version","dataset_id","registry_version","method_version","engine","correlation_engine","r_version","bioconductor_version","repositories","environment_lock","packages","session_info","bioc_valid","normalization","filter","annotation","composition","quantification","geoquery","contrasts","seeds","deg_rule","test_families","native_results","input_hashes","output_hashes"}
    for f in sorted(required-set(m)):report.error("METHOD_MISSING_FIELD",f,path)
    for f in required & set(m):
        if m[f] in (None,"",[],{}):report.error("METHOD_EMPTY_FIELD",f,path)
    if m.get("dataset_id")!=dataset_id:report.error("METHOD_DATASET_ID_MISMATCH","",path)
    if m.get("bioc_valid") is not True:report.error("BIOC_ENV_INVALID","BiocManager::valid must pass",path)
    env=m.get("environment_lock",{})
    if isinstance(env,dict) and env.get("type")=="renv":
        lock=safe_child(path.parent,str(env.get("path","")))
        if not lock or not lock.is_file() or sha(lock)!=env.get("sha256") or not valid_hash(env.get("sha256","")):report.error("ENVIRONMENT_LOCK_MISMATCH","",path)
        else:
            payload=load_json(lock,report) or {}
            if not all(k in payload for k in ("R","Repositories","Packages")):report.error("INVALID_RENV_LOCK","Missing R/Repositories/Packages",lock)
            else:
                if str(payload.get("R",{}).get("Version"))!=str(m.get("r_version")):report.error("LOCK_R_VERSION_MISMATCH","",lock)
                if str(payload.get("Bioconductor",{}).get("Version"))!=str(m.get("bioconductor_version")):report.error("LOCK_BIOC_VERSION_MISMATCH","",lock)
                lock_repos={x.get("Name"):x.get("URL") for x in payload.get("Repositories",[]) if isinstance(x,dict)}
                if lock_repos!=m.get("repositories"):report.error("LOCK_REPOSITORY_MISMATCH","",lock)
                lock_packages=payload.get("Packages",{})
                method_packages=m.get("packages",{})
                if not isinstance(lock_packages,dict):report.error("INVALID_RENV_LOCK","Packages must be an object",lock)
                if not isinstance(method_packages,dict):report.error("INVALID_METHOD_PACKAGES","packages must be an object",path);method_packages={}
                for pkg,version in method_packages.items():
                    if str(lock_packages.get(pkg,{}).get("Version"))!=str(version):report.error("LOCK_PACKAGE_MISMATCH",pkg,lock)
                base={"R","base","methods","stats","utils","graphics","grDevices","datasets","tools"}
                for pkg,record in lock_packages.items():
                    if not isinstance(record,dict) or any(missing(record.get(f)) for f in ("Package","Version","Source")) or record.get("Package")!=pkg:report.error("INVALID_RENV_PACKAGE_RECORD",pkg,lock);continue
                    for dep in renv_dependencies(record):
                        if dep not in base and dep not in lock_packages:report.error("LOCK_DEPENDENCY_MISSING",f"{pkg}->{dep}",lock)
    elif isinstance(env,dict) and env.get("type")=="container":
        if not re.fullmatch(r"sha256:[0-9a-f]{64}",str(env.get("digest",""))):report.error("INVALID_CONTAINER_DIGEST","",path)
    else:report.error("INVALID_ENVIRONMENT_LOCK_TYPE","",path)
    if isinstance(m.get("packages"),dict) and any(missing(v) for v in m["packages"].values()):report.error("EMPTY_PACKAGE_VERSION","",path)
    for section in ("input_hashes","output_hashes"):
        if isinstance(m.get(section),dict):
            for rel,v in m[section].items():
                p=safe_child(path.parent,rel)
                if not p or not p.is_file() or not valid_hash(v) or sha(p)!=v:report.error("METHOD_FILE_HASH_MISMATCH",f"{section}.{rel}",path)
    q=m.get("quantification",{})
    if isinstance(q,dict):
        needed={"source_type","quantifier","tx2gene_sha256","counts_from_abundance","offset_decision","constructor"}
        for f in needed-set(q):report.error("QUANTIFICATION_MISSING_FIELD",f,path)
        source=str(q.get("source_type","")).upper(); quantifier=q.get("quantifier"); cfa=q.get("counts_from_abundance"); offset=q.get("offset_decision"); constructor=q.get("constructor")
        engine=str(m.get("engine","")); combo=(source,quantifier,engine,constructor,cfa,offset)
        allowed={
            ("INTEGER_COUNTS","gene_counts","edgeR","DGEList","not_applicable","not_applicable"),
            ("INTEGER_COUNTS","gene_counts","DESeq2","DESeqDataSetFromMatrix","not_applicable","not_applicable"),
            ("INTEGER_COUNTS","gene_counts","limma_voom","DGEList_voom","not_applicable","not_applicable"),
            ("SALMON","salmon","DESeq2","DESeqDataSetFromTximport","no","average_transcript_length_preserved"),
            ("SALMON","salmon","edgeR","DGEList_with_offset","no","average_transcript_length_preserved"),
            ("SALMON","salmon","limma_voom","DGEList_voom","lengthScaledTPM","no_additional_length_offset"),
            ("RSEM","rsem","DESeq2","DESeqDataSetFromTximport","no","average_transcript_length_preserved"),
            ("RSEM","rsem","edgeR","DGEList_with_offset","no","average_transcript_length_preserved"),
            ("TPM_ONLY","expression_matrix","limma_continuous","ExpressionMatrix","not_applicable","not_applicable"),
            ("NORMALIZED_CONTINUOUS","expression_matrix","limma_continuous","ExpressionMatrix","not_applicable","not_applicable"),
            ("MICROARRAY_LOG2","array","limma","ExpressionSet","not_applicable","not_applicable"),
        }
        if combo not in allowed:report.error("INVALID_QUANTIFICATION_COMBINATION",repr(combo),path)
        if not valid_hash(q.get("tx2gene_sha256","")):report.error("INVALID_TX2GENE_HASH","",path)
    geo=m.get("geoquery",{})
    if isinstance(geo,dict):
        for f in ("package_version","calls","cache_manifest_sha256","platforms","alignment_verified","supplement_classification"):
            if f not in geo or geo[f] in (None,"",[],{}):report.error("GEOQUERY_MISSING_FIELD",f,path)
        if not valid_hash(geo.get("cache_manifest_sha256","")):report.error("INVALID_GEOQUERY_CACHE_HASH","",path)
        if geo.get("alignment_verified") is not True:report.error("GEOQUERY_ALIGNMENT_UNVERIFIED","",path)
    native_types=set();native_rows={}
    if isinstance(m.get("native_results"),list):
        for entry in m["native_results"]:
            if not isinstance(entry,dict):report.error("INVALID_NATIVE_RESULT_ENTRY","",path);continue
            p=safe_child(path.parent,str(entry.get("path","")));h=entry.get("sha256","");st=entry.get("statistic_type","")
            if not st or not p or not p.is_file() or not valid_hash(h) or sha(p)!=h:report.error("NATIVE_RESULT_HASH_MISMATCH",str(st),path)
            else:
                native_types.add(st); nr,_=read_tsv(p,report,{"dataset_id","model_id","contrast_id","estimand","test_id","feature_id","log2fc","statistic","df","p_value"})
                by={}
                for rn,row in enumerate(nr,2):
                    if row.get("dataset_id")!=dataset_id:report.error("NATIVE_DATASET_ID_MISMATCH","",p,rn)
                    tid=(row.get("model_id",""),row.get("contrast_id",""),row.get("estimand",""),row.get("test_id",""))
                    if tid in by:report.error("DUPLICATE_NATIVE_TEST_ID",repr(tid),p,rn)
                    by[tid]=row
                native_rows[st]=by
    families={}
    if isinstance(m.get("test_families"),list):
        for item in m["test_families"]:
            if not isinstance(item,dict) or not {"bh_family_id","model_id","contrast_id","estimand","test_universe_sha256"}.issubset(item):report.error("INVALID_METHOD_FAMILY","",path);continue
            fid=item["bh_family_id"]
            if fid in families:report.error("DUPLICATE_METHOD_FAMILY",fid,path)
            families[fid]=item
    m["_families"]=families
    m["_validated_native_types"]=native_types
    m["_native_rows"]=native_rows
    return m

def validate_de(path,report,dataset_id,meta,design,method):
    rows,cols=read_tsv(path,report,DE_REQ); seen=set(); feature_seen=set(); fam=defaultdict(list); expected_native=defaultdict(set)
    deg_rule=method.get("deg_rule",{}) if isinstance(method,dict) else {}; alpha=float(deg_rule.get("fdr_lte",0.05)); fc=float(deg_rule.get("abs_log2fc_gte",0.0))
    for rn,r in enumerate(rows,2):
        if r.get("dataset_id","").strip()!=dataset_id:report.error("DE_DATASET_ID_MISMATCH",f"Expected {dataset_id}",path,rn)
        test=r.get("test_id","").strip(); key=(r.get("model_id","").strip(),r.get("contrast_id","").strip(),r.get("estimand","").strip())
        test_key=key+(test,)
        if not test or test_key in seen:report.error("DUPLICATE_TEST_ID",test,path,rn)
        seen.add(test_key)
        feature_key=key+(r.get("feature_id","").strip(),)
        if not feature_key[-1] or feature_key in feature_seen:report.error("DUPLICATE_FEATURE_IN_FAMILY",repr(feature_key),path,rn)
        feature_seen.add(feature_key)
        if key not in design:report.error("UNKNOWN_DESIGN_FOREIGN_KEY",str(key),path,rn)
        fid=r.get("bh_family_id","");family=method.get("_families",{}).get(fid)
        if not family or (family.get("model_id"),family.get("contrast_id"),family.get("estimand"))!=key:report.error("METHOD_FAMILY_FOREIGN_KEY_MISMATCH",fid,path,rn)
        if r.get("annotation_version")!=str(method.get("annotation",{}).get("version")):report.error("ANNOTATION_VERSION_MISMATCH","",path,rn)
        if r.get("method_version")!=str(method.get("method_version")):report.error("METHOD_VERSION_MISMATCH","",path,rn)
        if family and r.get("test_universe_sha256")!=family.get("test_universe_sha256"):report.error("DE_UNIVERSE_FOREIGN_KEY_MISMATCH",fid,path,rn)
        if r.get("mapping_status") not in {"one_to_one","one_to_many","unmapped","retired","probe_collapsed"}:report.error("INVALID_MAPPING_STATUS",r.get("mapping_status",""),path,rn)
        status=r.get("test_status","")
        if status not in {"TESTED","FILTERED","COOKS_OUTLIER","MODEL_FAILED"}:report.error("INVALID_TEST_STATUS",status,path,rn);continue
        if status!="TESTED":
            if missing(r.get("not_tested_reason")):report.error("MISSING_NOT_TESTED_REASON",test,path,rn)
            if any(not missing(r.get(f)) for f in ("se","statistic","df","p_value","fdr")):report.error("UNTESTED_HAS_INFERENCE",test,path,rn)
            continue
        if not missing(r.get("not_tested_reason")):report.error("TESTED_HAS_NOT_TESTED_REASON",test,path,rn)
        vals={f:number(r.get(f,""),f,path,rn,report) for f in ("log2fc","se","statistic","p_value","fdr","base_expression")}
        if any(v is None for v in vals.values()):continue
        if vals["se"]<=0:report.error("INVALID_SE",test,path,rn)
        if not 0<=vals["p_value"]<=1 or not 0<=vals["fdr"]<=1:report.error("PROBABILITY_RANGE",test,path,rn)
        direction="up" if vals["log2fc"]>0 else "down" if vals["log2fc"]<0 else "zero"
        if r.get("direction")!=direction:report.error("DIRECTION_MISMATCH",test,path,rn)
        st=r.get("statistic_type","")
        df=number(r.get("df",""),"df",path,rn,report,allow_na=st=="wald_z")
        if st in {"moderated_t","t","ql_f"} and (df is None or df<=0):report.error("INVALID_STATISTIC_DF",test,path,rn)
        if st not in method.get("_validated_native_types",set()):report.error("UNVERIFIED_STATISTIC_TYPE",st,path,rn)
        expected_native[st].add(key+(test,))
        if st in {"wald_z","moderated_t","t"} and abs(vals["statistic"]-vals["log2fc"]/vals["se"])>1e-5:report.error("STATISTIC_MISMATCH",test,path,rn)
        if st=="wald_z" and abs(vals["p_value"]-math.erfc(abs(vals["statistic"])/math.sqrt(2)))>1e-6:report.error("P_VALUE_MISMATCH",test,path,rn)
        native=method.get("_native_rows",{}).get(st,{}).get(key+(test,))
        if not native:report.error("MISSING_NATIVE_TEST_ROW",test,path,rn)
        else:
            for f,observed in (("log2fc",vals["log2fc"]),("statistic",vals["statistic"]),("p_value",vals["p_value"])):
                nv=number(native.get(f,""),f,path,rn,report)
                if nv is not None and abs(nv-observed)>1e-8:report.error("NATIVE_RESULT_VALUE_MISMATCH",f"{test}:{f}",path,rn)
            if native.get("feature_id")!=r.get("feature_id"):report.error("NATIVE_FEATURE_ID_MISMATCH",test,path,rn)
            native_df=number(native.get("df",""),"df",path,rn,report,allow_na=st=="wald_z")
            if df is not None and native_df is not None and abs(df-native_df)>1e-8:report.error("NATIVE_RESULT_VALUE_MISMATCH",f"{test}:df",path,rn)
        observed=boolean(r.get("deg",""),"deg",path,rn,report); expected=vals["fdr"]<=alpha and abs(vals["log2fc"])>=fc
        if observed is not None and observed!=expected:report.error("DEG_RULE_MISMATCH",test,path,rn)
        if integer(r.get("n_samples",""),"n_samples",path,rn,report)!=meta["n_samples"]:report.error("COUNT_MISMATCH","n_samples",path,rn)
        if integer(r.get("n_donors",""),"n_donors",path,rn,report)!=meta["n_donors"]:report.error("COUNT_MISMATCH","n_donors",path,rn)
        if not valid_hash(r.get("test_universe_sha256","")):report.error("INVALID_TEST_UNIVERSE_HASH",test,path,rn)
        fam[r.get("bh_family_id","")].append((rn,vals["p_value"],vals["fdr"],r))
    for family,items in fam.items():
        if not family:report.error("EMPTY_BH_FAMILY","",path);continue
        for (rn,_,obs,_),exp in zip(items,bh([x[1] for x in items])):
            if abs(obs-exp)>max(1e-8,1e-6*exp):report.error("FDR_MISMATCH",family,path,rn)
    for st,native_rows in method.get("_native_rows",{}).items():
        if set(native_rows)!=expected_native.get(st,set()):report.error("NATIVE_RESULT_TEST_SET_MISMATCH",st,path)
    report.metrics[f"{dataset_id}_test_universe_rows"]=len(rows)
    all_by_family=defaultdict(list)
    for r in rows:all_by_family[r.get("bh_family_id","")].append(r)
    return {"tested":{family:[item[3].get("test_id","") for item in items] for family,items in fam.items()},"all":all_by_family}

def validate_source(path,report,dataset_id):
    rows,cols=read_tsv(path,report,SOURCE_REQ)
    for rn,r in enumerate(rows,2):
        if r.get("dataset_id")!=dataset_id:report.error("DATASET_ID_MISMATCH",dataset_id,path,rn)
        if not valid_hash(r.get("sha256","")):report.error("INVALID_SOURCE_HASH","Nonzero SHA-256 required",path,rn)
        if r.get("source_kind") not in {"GEO","SRA","BIOSAMPLE","PAPER","SUPPLEMENT","ANNOTATION","OTHER"}:report.error("INVALID_SOURCE_KIND","",path,rn)
        if r.get("content_role") not in {"counts","metadata","annotation","normalized_expression","transcript_estimates","raw_archive","method","other"}:report.error("INVALID_CONTENT_ROLE","",path,rn)
        if not valid_utc(r.get("retrieved_at_utc")) or not valid_utc(r.get("record_status_at_utc")):report.error("INVALID_UTC_TIMESTAMP","",path,rn)
        size=integer(r.get("bytes",""),"bytes",path,rn,report)
        if size is not None and size<0:report.error("NEGATIVE_SOURCE_SIZE","",path,rn)
        if r.get("access_status") not in {"PUBLIC","WITHDRAWN","RESTRICTED","UNKNOWN"}:report.error("INVALID_ACCESS_STATUS","",path,rn)
        if r.get("archive_check") not in {"PASS","NOT_APPLICABLE"}:report.error("INVALID_ARCHIVE_STATUS","",path,rn)
        if r.get("local_verification_status") not in {"VERIFIED","REMOTE_ONLY","UNAVAILABLE"}:report.error("INVALID_LOCAL_STATUS","",path,rn)
        if any(missing(r.get(f)) for f in ("record_status_at_utc","citation","terms_url","license_status","use_restrictions","privacy_handling")):report.error("INCOMPLETE_GOVERNANCE","",path,rn)
        exceptional=r.get("access_status")!="PUBLIC" or r.get("local_verification_status")!="VERIFIED"
        approved=(r.get("audit_exception_status")=="APPROVED" and all(not missing(r.get(f)) for f in ("audit_exception_id","audit_exception_approved_by","verification_reason","remote_immutable_id")))
        if exceptional and not approved:report.error("SOURCE_AUDIT_APPROVAL_REQUIRED",f"{r.get('access_status')}/{r.get('local_verification_status')}",path,rn)
        if r.get("access_status") in {"WITHDRAWN","RESTRICTED"}:report.error("SOURCE_NOT_PUBLISHABLE",r.get("access_status"),path,rn)
        if r.get("local_verification_status")=="VERIFIED":
            p=safe_child(path.parent,r.get("local_path",""))
            if not p or not p.is_file():report.error("LOCAL_SOURCE_MISSING","",path,rn)
            elif size!=p.stat().st_size or sha(p)!=r.get("sha256"):report.error("LOCAL_SOURCE_MISMATCH","",path,rn)

def validate_qc(path,report,dataset_id,expected_values):
    rows,_=read_tsv(path,report,QC_REQ)
    required={"included_samples","included_donors","excluded_samples","exclusion_reasons_complete","test_universe_sha256","design_rank","residual_df","normalization","mapping_coverage"}|set(expected_values); seen=set()
    for rn,r in enumerate(rows,2):
        if r.get("dataset_id")!=dataset_id:report.error("QC_DATASET_ID_MISMATCH","",path,rn)
        seen.add(r.get("metric",""))
        if r.get("status") not in {"PASS","FAIL","NOT_APPLICABLE"}:report.error("INVALID_QC_STATUS","",path,rn)
        if r.get("status")=="FAIL":report.error("QC_GATE_FAILED",r.get("metric",""),path,rn)
        metric=r.get("metric","")
        if metric in expected_values:
            expected=expected_values[metric]
            if isinstance(expected,(int,float)):
                value=number(r.get("value",""),"value",path,rn,report)
                if value is not None and abs(value-expected)>1e-8:report.error("QC_VALUE_MISMATCH",metric,path,rn)
            elif r.get("value")!=expected:report.error("QC_VALUE_MISMATCH",metric,path,rn)
    for m in required-seen:report.error("MISSING_QC_METRIC",m,path)

def validate_dataset_manifest(dataset_dir,dataset_id,report):
    path=dataset_dir/"dataset_manifest.json";m=load_json(path,report) or {}
    required={"schema_version","dataset_id","analysis_role","registry_version","canonical_directory","artifact_hashes","dataset_payload_sha256"}
    for f in required-set(m):report.error("DATASET_MANIFEST_MISSING_FIELD",f,path)
    expected_role=EXPECTED_DATASETS.get(dataset_id,{}).get("role")
    if m.get("dataset_id")!=dataset_id or m.get("canonical_directory")!=dataset_id or dataset_dir.name!=dataset_id:report.error("DATASET_MANIFEST_IDENTITY_MISMATCH",dataset_id,path)
    if m.get("analysis_role")!=expected_role:report.error("DATASET_MANIFEST_ROLE_MISMATCH",dataset_id,path)
    hashes=m.get("artifact_hashes",{}) if isinstance(m.get("artifact_hashes"),dict) else {}
    actual={p.relative_to(dataset_dir).as_posix() for p in dataset_dir.rglob("*") if p.is_file() and p.name!="dataset_manifest.json"}
    if set(hashes)!=actual:report.error("DATASET_MANIFEST_FILE_SET_MISMATCH",f"manifest={sorted(hashes)} actual={sorted(actual)}",path)
    for rel,h in hashes.items():
        p=safe_child(dataset_dir,rel)
        if not p or not p.is_file() or not valid_hash(h) or sha(p)!=h:report.error("DATASET_ARTIFACT_HASH_MISMATCH",rel,path)
    payload=hashlib.sha256("".join(f"{rel}\t{hashes[rel]}\n" for rel in sorted(hashes)).encode()).hexdigest()
    if m.get("dataset_payload_sha256")!=payload:report.error("DATASET_PAYLOAD_HASH_MISMATCH","",path)
    return m,sha(path) if path.is_file() else "",payload

def validate_dataset(dataset_dir,dataset_id,report):
    start_errors=report.error_count
    for f in ("dataset_manifest.json","sample_metadata.tsv","design.tsv","differential_expression.tsv","qc_summary.tsv","source_manifest.tsv","method.json","test_universe.tsv","feature_manifest.tsv"):
        if not (dataset_dir/f).is_file():report.error("MISSING_DATASET_ARTIFACT",f,dataset_dir/f)
    manifest,manifest_hash,dataset_payload=validate_dataset_manifest(dataset_dir,dataset_id,report)
    meta=validate_metadata(dataset_dir/"sample_metadata.tsv",report,dataset_id)
    if dataset_id=="GSE196728":
        clocks=sorted({r.get("time_of_day","").strip() for _,r in meta["included"] if not missing(r.get("time_of_day"))})
        expected={(donor,period,clock) for donor in meta["donors"] for period in ("sea_level","3800m","5100m") for clock in clocks}
        actual={(r.get("donor_id","").strip(),r.get("altitude_label","").strip(),r.get("time_of_day","").strip()) for _,r in meta["included"]}
        missing_cells=expected-actual; imbalance_path=dataset_dir/"planned_imbalance.tsv"
        declared=set()
        if missing_cells or imbalance_path.exists():
            ir,_=read_tsv(imbalance_path,report,{"planned_imbalance_id","donor_id","period","time_of_day","reason","model_handling","approved","approved_by"})
            record_ids=set()
            for rn,r in enumerate(ir,2):
                if any(missing(r.get(f)) for f in ("planned_imbalance_id","donor_id","period","time_of_day","reason","model_handling","approved_by")):report.error("INCOMPLETE_IMBALANCE_RECORD","",imbalance_path,rn)
                if boolean(r.get("approved",""),"approved",imbalance_path,rn,report) is not True:report.error("UNAPPROVED_IMBALANCE","",imbalance_path,rn)
                record_ids.add(r.get("planned_imbalance_id",""))
                declared.add((r.get("donor_id",""),r.get("period",""),r.get("time_of_day","")))
            if declared!=missing_cells:report.error("IMBALANCE_CELL_MISMATCH",f"declared={sorted(declared)} missing={sorted(missing_cells)}",imbalance_path)
            metadata_ids={r.get("planned_imbalance_id","").strip() for _,r in meta["included"] if not missing(r.get("planned_imbalance_id"))}
            if metadata_ids!=record_ids:report.error("IMBALANCE_ID_FOREIGN_KEY_MISMATCH","",imbalance_path)
    method=validate_method(dataset_dir/"method.json",report,dataset_id)
    # Dataset-specific sample flow is mandatory and blocks selective attrition.
    sf=method.get("sample_flow",{}) if isinstance(method,dict) else {}
    for f in ("expected_samples","expected_donors","minimum_donors","maximum_attrition","exclusion_reasons_complete","minimum_gate_rationale_path","minimum_gate_rationale_sha256"):
        if f not in sf:report.error("SAMPLE_FLOW_MISSING_FIELD",f,dataset_dir/"method.json")
    if sf:
        expected=EXPECTED_FLOW.get(dataset_id)
        if expected and (int(sf.get("expected_samples",-1)),int(sf.get("expected_donors",-1)))!=expected:report.error("EXPECTED_FLOW_MISMATCH",f"Expected original flow {expected}",dataset_dir/"method.json")
        if meta["n_donors"]<int(sf.get("minimum_donors",10**9)):report.error("BELOW_MINIMUM_DONORS","",dataset_dir)
        if int(sf.get("expected_samples",0))-meta["n_samples"]>int(sf.get("maximum_attrition",-1)):report.error("ATTRITION_EXCEEDS_LIMIT","",dataset_dir)
        if sf.get("exclusion_reasons_complete") is not True:report.error("INCOMPLETE_EXCLUSION_FLOW","",dataset_dir)
        rationale=safe_child(dataset_dir,str(sf.get("minimum_gate_rationale_path","")))
        if not rationale or not rationale.is_file() or sha(rationale)!=sf.get("minimum_gate_rationale_sha256") or not valid_hash(sf.get("minimum_gate_rationale_sha256","")):report.error("MINIMUM_GATE_RATIONALE_MISMATCH","",dataset_dir/"method.json")
    design=validate_design(dataset_dir/"design.tsv",report,dataset_id,meta)
    for family_id,family in method.get("_families",{}).items():
        family_key=(family.get("model_id",""),family.get("contrast_id",""),family.get("estimand",""))
        if family_key not in design:report.error("METHOD_FAMILY_DESIGN_FOREIGN_KEY_MISMATCH",family_id,dataset_dir/"method.json")
    correlation_engines={r.get("correlation_engine") for r in design.values()}
    if len(correlation_engines)!=1 or method.get("correlation_engine") not in correlation_engines:report.error("METHOD_CORRELATION_ENGINE_MISMATCH","",dataset_dir/"method.json")
    de_result=validate_de(dataset_dir/"differential_expression.tsv",report,dataset_id,meta,design,method)
    validate_source(dataset_dir/"source_manifest.tsv",report,dataset_id)
    features,_=read_tsv(dataset_dir/"feature_manifest.tsv",report,{"dataset_id","feature_id","feature_namespace","feature_namespace_version"})
    feature_ids=[]
    for rn,r in enumerate(features,2):
        if r.get("dataset_id")!=dataset_id:report.error("FEATURE_DATASET_ID_MISMATCH","",dataset_dir/"feature_manifest.tsv",rn)
        if not r.get("feature_id") or r.get("feature_id") in feature_ids:report.error("DUPLICATE_INPUT_FEATURE",r.get("feature_id",""),dataset_dir/"feature_manifest.tsv",rn)
        feature_ids.append(r.get("feature_id",""))
    feature_hash=sha(dataset_dir/"feature_manifest.tsv") if (dataset_dir/"feature_manifest.tsv").is_file() else ""
    universe,_=read_tsv(dataset_dir/"test_universe.tsv",report,{"dataset_id","model_id","contrast_id","estimand","bh_family_id","input_feature_count","tested_count","input_feature_sha256","tested_feature_sha256","family_universe_sha256"})
    universe_families=set()
    for rn,r in enumerate(universe,2):
        family=r.get("bh_family_id","");universe_families.add(family)
        if r.get("dataset_id")!=dataset_id:report.error("UNIVERSE_DATASET_ID_MISMATCH","",dataset_dir/"test_universe.tsv",rn)
        mf=method.get("_families",{}).get(family,{})
        if (r.get("model_id"),r.get("contrast_id"),r.get("estimand"))!=(mf.get("model_id"),mf.get("contrast_id"),mf.get("estimand")):report.error("UNIVERSE_FAMILY_FOREIGN_KEY_MISMATCH",family,dataset_dir/"test_universe.tsv",rn)
        for f in ("input_feature_sha256","tested_feature_sha256"):
            if not valid_hash(r.get(f,"")):report.error("INVALID_UNIVERSE_HASH",f,dataset_dir/"test_universe.tsv",rn)
        input_count=integer(r.get("input_feature_count",""),"input_feature_count",dataset_dir/"test_universe.tsv",rn,report)
        if input_count!=len(feature_ids) or r.get("input_feature_sha256")!=feature_hash:report.error("INPUT_FEATURE_UNIVERSE_MISMATCH",family,dataset_dir/"test_universe.tsv",rn)
        count=integer(r.get("tested_count",""),"tested_count",dataset_dir/"test_universe.tsv",rn,report)
        ids=de_result.get("tested",{}).get(family,[])
        expected_hash=hashlib.sha256(("\n".join(ids)+"\n").encode()).hexdigest()
        if count!=len(ids):report.error("TESTED_COUNT_MISMATCH",family,dataset_dir/"test_universe.tsv",rn)
        if r.get("tested_feature_sha256")!=expected_hash:report.error("TESTED_UNIVERSE_HASH_MISMATCH",family,dataset_dir/"test_universe.tsv",rn)
        family_rows=de_result.get("all",{}).get(family,[])
        by_feature={x.get("feature_id"):x for x in family_rows}
        if set(by_feature)!=set(feature_ids) or len(family_rows)!=len(feature_ids):report.error("INCOMPLETE_FEATURE_UNIVERSE",family,dataset_dir/"differential_expression.tsv")
        family_payload="".join(f"{fid}\t{by_feature.get(fid,{}).get('test_id','')}\t{by_feature.get(fid,{}).get('test_status','')}\n" for fid in feature_ids)
        family_hash=hashlib.sha256(family_payload.encode()).hexdigest()
        if r.get("family_universe_sha256")!=family_hash or mf.get("test_universe_sha256")!=family_hash:report.error("FAMILY_UNIVERSE_HASH_MISMATCH",family,dataset_dir/"test_universe.tsv",rn)
        for x in family_rows:
            if x.get("test_universe_sha256")!=family_hash:report.error("DE_UNIVERSE_FOREIGN_KEY_MISMATCH",family,dataset_dir/"differential_expression.tsv")
    for family in set(de_result.get("all",{}))-universe_families:report.error("MISSING_TEST_UNIVERSE_FAMILY",family,dataset_dir/"test_universe.tsv")
    expected_qc={"included_samples":meta["n_samples"],"included_donors":meta["n_donors"],"excluded_samples":len(meta["rows"])-meta["n_samples"],"input_feature_count":len(feature_ids),"tested_count_total":sum(len(x) for x in de_result.get("tested",{}).values()),"design_rank_min":min((int(x.get("design_rank",0)) for x in design.values()),default=0),"residual_df_min":min((int(x.get("residual_df",0)) for x in design.values()),default=0),"test_universe_sha256":sha(dataset_dir/"test_universe.tsv") if (dataset_dir/"test_universe.tsv").is_file() else ""}
    validate_qc(dataset_dir/"qc_summary.tsv",report,dataset_id,expected_qc)
    if dataset_id in {"GSE196728","GSE75665","GSE103940"}:
        estimands={k[2] for k in design}
        if not {"total_bulk","composition_conditional"}.issubset(estimands):report.error("MISSING_WHOLE_BLOOD_ESTIMAND","Both total_bulk and composition_conditional are required",dataset_dir/"design.tsv")
    return {"dataset_id":dataset_id,"valid":report.error_count==start_errors,"metadata":meta,"design":design,"method":method,"de":de_result,"manifest":manifest,"manifest_sha256":manifest_hash,"dataset_payload_sha256":dataset_payload,"feature_hash":feature_hash,"feature_ids":feature_ids,"universe_rows":universe}

def validate_registry(project_root,metadata,report):
    rows,cols=read_tsv(project_root/"config/expression_datasets.tsv",report,{"dataset_id","tissue","design_role","analysis_role","de_authorized"})
    by={}
    for rn,r in enumerate(rows,2):
        ds=r.get("dataset_id","")
        if ds in by:report.error("DUPLICATE_REGISTRY_DATASET_ID",ds,project_root/"config/expression_datasets.tsv",rn)
        by[ds]=r
    for ds,exp in EXPECTED_DATASETS.items():
        if ds not in by:report.error("ACTIVE_DATASET_MISSING",ds,project_root/"config/expression_datasets.tsv");continue
        if by[ds].get("analysis_role")!=exp["role"]:report.error("REGISTRY_ROLE_MISMATCH",ds,project_root/"config/expression_datasets.tsv")
    for ds in EXPECTED_DATASETS:
        if ds in by and boolean(by[ds].get("de_authorized",""),"de_authorized",project_root/"config/expression_datasets.tsv",0,report) is not True:
            report.error("ACTIVE_DE_NOT_AUTHORIZED",ds,project_root/"config/expression_datasets.tsv")
    if metadata:validate_metadata(metadata,report)

def find_datasets(root,report):
    out={};manifest_dirs=set()
    for p in root.rglob("dataset_manifest.json"):
        manifest_dirs.add(p.parent.resolve());m=load_json(p,report) or {};ds=str(m.get("dataset_id","")).strip()
        if ds not in AUTHORIZED_DATASETS:report.error("UNAUTHORIZED_DATASET_ID",ds,p)
        if ds in out:report.error("DUPLICATE_DATASET_ID",ds,p)
        else:out[ds]=p.parent
    scoped={"sample_metadata.tsv","design.tsv","differential_expression.tsv","method.json","qc_summary.tsv","source_manifest.tsv","test_universe.tsv","feature_manifest.tsv","cell_composition.tsv","planned_imbalance.tsv","design_matrix.tsv","contrast.tsv","block_vector.tsv","native_results.tsv"}
    for p in root.rglob("*"):
        if p.is_file() and p.name in scoped and p.parent.resolve() not in manifest_dirs:report.error("ORPHAN_DATASET_ARTIFACT",p.name,p)
    return out
def validate_eligibility(path,report):
    rows,_=read_tsv(path,report,ELIG_REQ); by={}
    for rn,r in enumerate(rows,2):
        ds=r.get("dataset_id","")
        if ds not in AUTHORIZED_DATASETS:report.error("UNAUTHORIZED_ELIGIBILITY_DATASET",ds,path,rn)
        if ds in by:report.error("DUPLICATE_ELIGIBILITY_DATASET",ds,path,rn)
        by[ds]=r
        if r.get("eligibility_status") not in {"ELIGIBLE","INELIGIBLE","CONDITIONAL","NOT_TRIGGERED"}:report.error("INVALID_ELIGIBILITY_STATUS","",path,rn)
        role=r.get("analysis_role",""); auth=boolean(r.get("de_authorized",""),"de_authorized",path,rn,report)
        gates=("pairing_complete","donor_independence_verified","tissue_match","primary_contrast_estimable","design_full_rank","residual_df_adequate","source_data_available","exclusion_reasons_complete")
        parsed={g:boolean(r.get(g,""),g,path,rn,report) for g in gates}
        if r.get("eligibility_status")=="ELIGIBLE" and (auth is not True or any(v is not True for v in parsed.values())):report.error("FALSE_ELIGIBLE_CLAIM",ds,path,rn)
        n=integer(r.get("n_donors_included",""),"n_donors_included",path,rn,report); minimum=integer(r.get("minimum_donors",""),"minimum_donors",path,rn,report)
        expected=EXPECTED_FLOW.get(ds)
        if expected:
            es=integer(r.get("expected_samples",""),"expected_samples",path,rn,report);ed=integer(r.get("expected_donors",""),"expected_donors",path,rn,report)
            if (es,ed)!=expected:report.error("EXPECTED_FLOW_MISMATCH",ds,path,rn)
        if r.get("eligibility_status")=="ELIGIBLE" and n is not None and minimum is not None and n<minimum:report.error("BELOW_MINIMUM_DONORS",ds,path,rn)
        expected_role=EXPECTED_DATASETS.get(ds,{}).get("role")
        if expected_role and role!=expected_role:report.error("ELIGIBILITY_ROLE_MISMATCH",ds,path,rn)
    for ds in EXPECTED_DATASETS:
        if ds not in by:report.error("MISSING_DATASET_ROW",ds,path)
    return by
def canonical_index(rows):return "".join(f"{r['path']}\t{r['bytes']}\t{r['sha256']}\n" for r in sorted(rows,key=lambda x:x["path"]))
def validate_artifact_index(root,report):
    path=root/"artifact_index.tsv"; rows,_=read_tsv(path,report,{"path","bytes","sha256"}); seen=set(); indexed=set()
    for rn,r in enumerate(rows,2):
        rel=r.get("path",""); p=safe_child(root,rel)
        if not p:report.error("UNSAFE_ARTIFACT_PATH",rel,path,rn);continue
        if rel in seen:report.error("DUPLICATE_ARTIFACT_PATH",rel,path,rn)
        seen.add(rel);indexed.add(rel)
        if not p.is_file() or p.is_symlink():report.error("INDEXED_ARTIFACT_INVALID",rel,path,rn);continue
        size=integer(r.get("bytes",""),"bytes",path,rn,report)
        if size!=p.stat().st_size or sha(p)!=r.get("sha256") or not valid_hash(r.get("sha256","")):report.error("ARTIFACT_HASH_MISMATCH",rel,path,rn)
    excluded={"artifact_index.tsv","validation_report.json","_SUCCESS.json"}
    actual={p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() and p.relative_to(root).as_posix() not in excluded}
    for rel in sorted(actual-indexed):report.error("UNINDEXED_ARTIFACT",rel,path)
    for rel in sorted(indexed-actual):report.error("STALE_ARTIFACT_INDEX",rel,path)
    index_hash=sha(path) if path.is_file() else ""; payload=hashlib.sha256(canonical_index(rows).encode()).hexdigest()
    return index_hash,payload

def close_number(observed,expected,tol=1e-8):
    return abs(observed-expected)<=max(tol,tol*abs(expected))

def composition_values(total,conditional):
    delta=conditional-total
    if abs(total)<1e-12:
        return delta,None,"stable_zero","unresolved"
    attenuation=1-abs(conditional)/abs(total)
    if total*conditional<0:return delta,attenuation,"direction_changed","direction_changed"
    relation="conditional_zero" if abs(conditional)<1e-12 else "same_direction"
    if attenuation>=0.5:return delta,attenuation,relation,"attenuated"
    if abs(delta)>=0.5:return delta,attenuation,relation,"composition_sensitive"
    return delta,attenuation,relation,"robust_residual"

def source_effects_sha256(rows):
    payload=""
    for r in sorted(rows,key=lambda x:x.get("dataset_id","")):
        effect=float(r["effect"]);se=float(r["se"])
        payload+=f"{r['dataset_id']}\t{r['source_family_id']}\t{r['model_id']}\t{r['contrast_id']}\t{r['scientific_comparison_id']}\t{r['estimand']}\t{r['test_id']}\t{effect:.17g}\t{se:.17g}\t{r['source_de_sha256']}\n"
    return hashlib.sha256(payload.encode()).hexdigest()

def reml_tau2(effects,ses):
    variances=[x*x for x in ses]
    def objective(tau2):
        weights=[1/(v+tau2) for v in variances];mu=sum(w*y for w,y in zip(weights,effects))/sum(weights)
        return sum(math.log(v+tau2) for v in variances)+math.log(sum(weights))+sum(w*(y-mu)**2 for w,y in zip(weights,effects))
    upper=max(max(variances),max(effects)-min(effects),1e-6)
    while objective(upper)<objective(upper/2) and upper<1e6:upper*=2
    lo,hi=0.0,upper;phi=(1+math.sqrt(5))/2
    for _ in range(100):
        c=hi-(hi-lo)/phi;d=lo+(hi-lo)/phi
        if objective(c)<objective(d):hi=d
        else:lo=c
    estimate=(lo+hi)/2
    return 0.0 if objective(0)<=objective(estimate)+1e-12 else estimate

def meta_statistics(source_rows):
    effects=[float(r["effect"]) for r in source_rows];ses=[float(r["se"]) for r in source_rows];k=len(effects);tau2=reml_tau2(effects,ses)
    weights=[1/(s*s+tau2) for s in ses];mu=sum(w*y for w,y in zip(weights,effects))/sum(weights)
    qscale=sum(w*(y-mu)**2 for w,y in zip(weights,effects))/(k-1);se_hksj=math.sqrt(qscale/sum(weights));se_modified=math.sqrt(max(1.0,qscale)/sum(weights))
    statistic=mu/se_hksj if se_hksj>0 else math.copysign(math.inf,mu) if mu else 0.0
    if k==2:p_value=1-2/math.pi*math.atan(abs(statistic))
    else:p_value=float("nan")
    fixed_weights=[1/(s*s) for s in ses];fixed=sum(w*y for w,y in zip(fixed_weights,effects))/sum(fixed_weights);fixed_se=math.sqrt(1/sum(fixed_weights));q=sum(w*(y-fixed)**2 for w,y in zip(fixed_weights,effects));i2=max(0.0,(q-(k-1))/q)*100 if q>0 else 0.0
    concordant=all(y>0 for y in effects) or all(y<0 for y in effects) or all(abs(y)<1e-12 for y in effects)
    tcrit=12.706204736432095 if k==2 else float("nan")
    return {"meta_log2fc":mu,"meta_se":se_hksj,"meta_statistic":statistic,"df":k-1,"p_value":p_value,"q_heterogeneity":q,"i2":i2,"tau2":tau2,"direction_concordant":"true" if concordant else "false","modified_hksj_low":mu-tcrit*se_modified,"modified_hksj_high":mu+tcrit*se_modified,"fixed_effect_log2fc":fixed,"fixed_effect_se":fixed_se}

def validate_composition(root,dirs,summaries,eligible,report):
    path=root/"cell_composition_sensitivity.tsv";rows,_=read_tsv(path,report,COMPOSITION_REQ)
    allowed={ds for ds,s in summaries.items() if s.get("manifest",{}).get("analysis_role")=="WHOLE_BLOOD_DE" and s.get("valid")}
    de_by_dataset={}
    for ds in allowed:
        de_rows,_=read_tsv(dirs[ds]/"differential_expression.tsv",Report("join"),DE_REQ)
        de_by_dataset[ds]={(r.get("model_id"),r.get("contrast_id"),r.get("estimand"),r.get("test_id")):r for r in de_rows if r.get("test_status")=="TESTED"}
    observed=defaultdict(set);seen=set()
    for rn,r in enumerate(rows,2):
        ds=r.get("dataset_id","");de_hash=sha(dirs[ds]/"differential_expression.tsv") if ds in dirs and (dirs[ds]/"differential_expression.tsv").is_file() else ""
        if ds not in allowed or r.get("source_de_sha256")!=de_hash:report.error("COMPOSITION_SOURCE_FOREIGN_KEY_MISMATCH",ds,path,rn);continue
        total_key=(r.get("total_model_id"),r.get("contrast_id"),"total_bulk",r.get("test_id"))
        conditional_key=(r.get("conditional_model_id"),r.get("contrast_id"),"composition_conditional",r.get("test_id"))
        row_key=(ds,total_key,conditional_key)
        if row_key in seen:report.error("DUPLICATE_COMPOSITION_TEST",repr(row_key),path,rn)
        seen.add(row_key)
        total=de_by_dataset[ds].get(total_key);conditional=de_by_dataset[ds].get(conditional_key)
        if not total or not conditional:report.error("COMPOSITION_DE_FOREIGN_KEY_MISMATCH",repr((total_key,conditional_key)),path,rn);continue
        observed[ds].add(r.get("test_id"))
        pairs=(("total_bulk_log2fc",total["log2fc"]),("total_bulk_se",total["se"]),("total_bulk_statistic",total["statistic"]),("total_bulk_fdr",total["fdr"]),("conditional_log2fc",conditional["log2fc"]),("conditional_se",conditional["se"]),("conditional_statistic",conditional["statistic"]),("conditional_fdr",conditional["fdr"]))
        parsed={}
        for field,expected_text in pairs:
            value=number(r.get(field,""),field,path,rn,report);expected=float(expected_text)
            parsed[field]=value
            if value is not None and not close_number(value,expected):report.error("COMPOSITION_VALUE_MISMATCH",field,path,rn)
        if parsed.get("total_bulk_log2fc") is None or parsed.get("conditional_log2fc") is None:continue
        delta,attenuation,relation,sensitivity=composition_values(float(total["log2fc"]),float(conditional["log2fc"]))
        observed_delta=number(r.get("delta_log2fc",""),"delta_log2fc",path,rn,report)
        observed_attenuation=number(r.get("attenuation_fraction",""),"attenuation_fraction",path,rn,report,allow_na=attenuation is None)
        if observed_delta is not None and not close_number(observed_delta,delta):report.error("COMPOSITION_DERIVATION_MISMATCH","delta_log2fc",path,rn)
        if (attenuation is None and observed_attenuation is not None) or (attenuation is not None and (observed_attenuation is None or not close_number(observed_attenuation,attenuation))):report.error("COMPOSITION_DERIVATION_MISMATCH","attenuation_fraction",path,rn)
        if r.get("direction_relation")!=relation or r.get("sensitivity_class")!=sensitivity:report.error("COMPOSITION_CLASS_MISMATCH",f"{relation}/{sensitivity}",path,rn)
    for ds in eligible:
        source_ids={k[3] for k in de_by_dataset.get(ds,{}) if k[2]=="total_bulk"} & {k[3] for k in de_by_dataset.get(ds,{}) if k[2]=="composition_conditional"}
        if observed.get(ds,set())!=source_ids:report.error("COMPOSITION_TEST_UNIVERSE_MISMATCH",ds,path)

def validate_pbmc(root,dirs,summaries,eligible,report):
    hypothesis_path=root/"whole_blood_hypotheses.tsv";hypotheses,_=read_tsv(hypothesis_path,report,HYPOTHESIS_REQ)
    hypothesis_hash=sha(hypothesis_path) if hypothesis_path.is_file() else "";by_target={}
    for rn,r in enumerate(hypotheses,2):
        target=r.get("target_id","")
        if any(missing(r.get(f)) for f in HYPOTHESIS_REQ):report.error("INCOMPLETE_HYPOTHESIS_TARGET",target,hypothesis_path,rn)
        if not target or target in by_target:report.error("DUPLICATE_HYPOTHESIS_TARGET",target,hypothesis_path,rn)
        by_target[target]=r;ds=r.get("source_dataset_id","")
        if ds not in eligible:report.error("INVALID_HYPOTHESIS_SOURCE_DATASET",ds,hypothesis_path,rn);continue
        if not valid_utc(r.get("frozen_at_utc")) or r.get("registry_version")!=summaries.get(ds,{}).get("manifest",{}).get("registry_version"):report.error("INVALID_HYPOTHESIS_FREEZE","",hypothesis_path,rn)
        de_path=dirs[ds]/"differential_expression.tsv"
        if r.get("source_de_sha256")!=sha(de_path):report.error("HYPOTHESIS_SOURCE_HASH_MISMATCH","",hypothesis_path,rn)
        de_rows,_=read_tsv(de_path,Report("hypothesis-join"),DE_REQ)
        source=next((x for x in de_rows if x.get("model_id")==r.get("source_model_id") and x.get("contrast_id")==r.get("source_contrast_id") and x.get("estimand")==r.get("source_estimand") and x.get("test_id")==r.get("source_test_id") and x.get("test_status")=="TESTED"),None)
        if not source:report.error("HYPOTHESIS_SOURCE_TEST_MISSING",target,hypothesis_path,rn)
    pbmc_path=root/"altitude_pbmc_validation.tsv";pbmc_rows,_=read_tsv(pbmc_path,report,PBMC_REQ)
    pbmc_de_path=dirs.get("GSE46480",Path("/nonexistent"))/"differential_expression.tsv"
    pbmc_de,_=read_tsv(pbmc_de_path,Report("pbmc-join"),DE_REQ)
    pbmc_by_key={(r.get("model_id"),r.get("contrast_id"),r.get("estimand"),r.get("test_id")):r for r in pbmc_de if r.get("test_status")=="TESTED"}
    seen=set()
    for rn,r in enumerate(pbmc_rows,2):
        target=r.get("target_id","");hyp=by_target.get(target);key=(r.get("model_id"),r.get("contrast_id"),r.get("estimand"),r.get("test_id"));source=pbmc_by_key.get(key)
        if target in seen:report.error("DUPLICATE_PBMC_TARGET",target,pbmc_path,rn)
        seen.add(target)
        if r.get("dataset_id")!="GSE46480" or r.get("hypothesis_status")!="PRESPECIFIED_TARGET" or not hyp:report.error("INVALID_PBMC_TARGET_ROW","",pbmc_path,rn);continue
        expected_key=(hyp.get("pbmc_model_id"),hyp.get("pbmc_contrast_id"),hyp.get("pbmc_estimand"),hyp.get("pbmc_test_id"))
        if key!=expected_key or r.get("multiplicity_family_id")!=hyp.get("multiplicity_family_id") or r.get("registry_version")!=hyp.get("registry_version") or r.get("target_source_sha256")!=hypothesis_hash:report.error("PBMC_TARGET_FOREIGN_KEY_MISMATCH",target,pbmc_path,rn)
        if not source or source.get("bh_family_id")!=r.get("multiplicity_family_id"):report.error("PBMC_DE_FOREIGN_KEY_MISMATCH",repr(key),pbmc_path,rn);continue
        for field in ("statistic","p_value","fdr"):
            value=number(r.get(field,""),field,pbmc_path,rn,report);expected=float(source[field])
            if value is not None and not close_number(value,expected):report.error("PBMC_VALUE_MISMATCH",field,pbmc_path,rn)
    if seen!=set(by_target):report.error("PBMC_TARGET_SET_MISMATCH","",pbmc_path)

def validate_meta(root,dirs,summaries,eligible,report):
    meta_path=root/"altitude_whole_blood_meta.tsv";mm_path=root/"meta_method.json";effects_path=root/"meta_study_effects.tsv"
    meta_rows,_=read_tsv(meta_path,report,META_RESULT_REQ);effects,_=read_tsv(effects_path,report,META_EFFECT_REQ);mm=load_json(mm_path,report) or {}
    required_meta={"k":2,"tau2_estimator":"REML","interval_method":"HKSJ","interpretation":"DESCRIPTIVE","confirmatory_deg":False}
    if any(mm.get(k)!=v for k,v in required_meta.items()) or not {"modified_HKSJ","fixed_effect"}.issubset(set(mm.get("sensitivity_methods",[]))):report.error("INVALID_K2_META_METHOD","",mm_path)
    contributors=mm.get("contributors",[]);contributor_ids=[x.get("dataset_id") for x in contributors if isinstance(x,dict)]
    if set(contributor_ids)!=set(eligible) or len(contributor_ids)!=2:report.error("META_CONTRIBUTOR_SET_MISMATCH",repr(contributor_ids),mm_path)
    selections={};source_de={}
    for entry in contributors:
        if not isinstance(entry,dict):report.error("INVALID_META_CONTRIBUTOR","",mm_path);continue
        required={"dataset_id","source_family_id","model_id","contrast_id","estimand","scientific_comparison_id","test_universe_sha256","source_de_sha256","dataset_payload_sha256"}
        if any(missing(entry.get(f)) for f in required):report.error("INCOMPLETE_META_SOURCE_SELECTION",repr(entry),mm_path);continue
        ds=entry["dataset_id"];summary=summaries.get(ds,{});de_path=dirs.get(ds,Path("/nonexistent"))/"differential_expression.tsv";family=summary.get("method",{}).get("_families",{}).get(entry["source_family_id"])
        design=summary.get("design",{}).get((entry["model_id"],entry["contrast_id"],entry["estimand"]))
        if not de_path.is_file() or entry["source_de_sha256"]!=sha(de_path) or entry["dataset_payload_sha256"]!=summary.get("dataset_payload_sha256"):report.error("META_CONTRIBUTOR_HASH_MISMATCH",ds,mm_path)
        if not family or (family.get("model_id"),family.get("contrast_id"),family.get("estimand"))!=(entry["model_id"],entry["contrast_id"],entry["estimand"]) or family.get("test_universe_sha256")!=entry["test_universe_sha256"]:report.error("META_SOURCE_FAMILY_MISMATCH",ds,mm_path)
        if not design or design.get("scientific_comparison_id")!=entry["scientific_comparison_id"] or entry["estimand"]!="total_bulk":report.error("META_SOURCE_COMPARISON_MISMATCH",ds,mm_path)
        selections[ds]=entry
        de_rows,_=read_tsv(de_path,Report("meta-join"),DE_REQ)
        source_de[ds]={(r.get("model_id"),r.get("contrast_id"),r.get("estimand"),r.get("test_id")):r for r in de_rows if r.get("bh_family_id")==entry["source_family_id"] and r.get("test_status")=="TESTED"}
    comparisons={x.get("scientific_comparison_id") for x in selections.values()}
    if len(comparisons)!=1:report.error("META_CROSS_STUDY_COMPARISON_MISMATCH",repr(comparisons),mm_path)
    meta_family=mm.get("meta_family_id");comparison=next(iter(comparisons),None)
    if missing(meta_family) or mm.get("scientific_comparison_id")!=comparison or mm.get("estimand")!="total_bulk":report.error("INVALID_META_FAMILY_DECLARATION","",mm_path)
    effect_by_key={};study_ids={}
    for rn,r in enumerate(effects,2):
        ds=r.get("dataset_id","");selection=selections.get(ds);full_key=(r.get("model_id"),r.get("contrast_id"),r.get("estimand"),r.get("test_id"))
        if not selection or any(r.get(f)!=selection.get(f) for f in ("source_family_id","model_id","contrast_id","scientific_comparison_id","estimand")):report.error("META_STUDY_EFFECT_FOREIGN_KEY_MISMATCH","",effects_path,rn);continue
        source=source_de.get(ds,{}).get(full_key)
        if not source or r.get("source_de_sha256")!=selection.get("source_de_sha256"):report.error("META_STUDY_EFFECT_FOREIGN_KEY_MISMATCH","",effects_path,rn);continue
        observed_effect=number(r.get("effect",""),"effect",effects_path,rn,report);observed_se=number(r.get("se",""),"se",effects_path,rn,report)
        if observed_effect is not None and observed_se is not None and (not close_number(observed_effect,float(source["log2fc"])) or not close_number(observed_se,float(source["se"]))):report.error("META_STUDY_EFFECT_VALUE_MISMATCH","",effects_path,rn)
        key=(ds,)+full_key
        if key in effect_by_key:report.error("DUPLICATE_META_STUDY_EFFECT",repr(key),effects_path,rn)
        effect_by_key[key]=r;study_ids.setdefault(ds,set()).add(r.get("test_id"))
    source_sets={ds:{k[3] for k in source_de.get(ds,{})} for ds in eligible}
    if any(study_ids.get(ds,set())!=source_sets.get(ds,set()) for ds in eligible) or len({frozenset(v) for v in source_sets.values()})!=1:report.error("META_STUDY_EFFECT_COVERAGE_MISMATCH",repr(study_ids),effects_path)
    expected_ids=set.intersection(*source_sets.values()) if source_sets else set()
    native_spec=mm.get("native_results",{});native_path=safe_child(root,str(native_spec.get("path",""))) if isinstance(native_spec,dict) else None
    if not native_path or not native_path.is_file() or sha(native_path)!=native_spec.get("sha256") or any(missing(native_spec.get(f)) for f in ("package","package_version","statistic_contract")) or native_spec.get("statistic_contract")!="REML_HKSJ":report.error("META_NATIVE_ARTIFACT_INVALID","",mm_path);native_rows=[]
    else:native_rows,_=read_tsv(native_path,report,META_RESULT_REQ)
    native_by={r.get("test_id"):r for r in native_rows};published_by={r.get("test_id"):r for r in meta_rows}
    if set(native_by)!=expected_ids or len(native_rows)!=len(expected_ids) or set(published_by)!=expected_ids or len(meta_rows)!=len(expected_ids):report.error("META_RESULT_TEST_SET_MISMATCH","",meta_path)
    numeric_fields=("meta_log2fc","meta_se","meta_statistic","df","p_value","fdr","q_heterogeneity","i2","tau2","modified_hksj_low","modified_hksj_high","fixed_effect_log2fc","fixed_effect_se")
    for test_id in expected_ids:
        source_rows=[effect_by_key[(ds,selections[ds]["model_id"],selections[ds]["contrast_id"],selections[ds]["estimand"],test_id)] for ds in eligible if (ds,selections[ds]["model_id"],selections[ds]["contrast_id"],selections[ds]["estimand"],test_id) in effect_by_key]
        source_features=[]
        for ds in eligible:
            source=source_de.get(ds,{}).get((selections[ds]["model_id"],selections[ds]["contrast_id"],selections[ds]["estimand"],test_id))
            if source:source_features.append((source.get("feature_id"),source.get("feature_namespace"),source.get("feature_namespace_version")))
        if len(set(source_features))!=1:report.error("META_SOURCE_FEATURE_MISMATCH",f"{test_id}: {source_features}",effects_path)
        native=native_by.get(test_id,{});published=published_by.get(test_id,{})
        expected_source_hash=source_effects_sha256(source_rows) if len(source_rows)==len(eligible) else ""
        identity={"meta_family_id":meta_family,"scientific_comparison_id":comparison,"estimand":"total_bulk","test_id":test_id,"k_studies":"2","estimand_comparability_status":"COMPARABLE","source_effects_sha256":expected_source_hash}
        for field,expected in identity.items():
            if native.get(field)!=expected or published.get(field)!=expected:report.error("META_RESULT_IDENTITY_MISMATCH",f"{test_id}:{field}",meta_path)
        for field in numeric_fields:
            nv=number(native.get(field,""),field,native_path or mm_path,0,report);pv=number(published.get(field,""),field,meta_path,0,report)
            if nv is not None and pv is not None and not close_number(pv,nv):report.error("META_RESULT_VALUE_MISMATCH",f"{test_id}:{field}",meta_path)
        if len(source_rows)==2 and all(float(x["se"])>0 for x in source_rows):
            derived=meta_statistics(source_rows)
            for field,expected in derived.items():
                if field=="direction_concordant":
                    if native.get(field)!=expected or published.get(field)!=expected:report.error("META_DIRECTION_MISMATCH",test_id,meta_path)
                else:
                    nv=number(native.get(field,""),field,native_path or mm_path,0,report)
                    if nv is not None and not close_number(nv,expected):report.error("META_NATIVE_DERIVATION_MISMATCH",f"{test_id}:{field}",native_path or mm_path)
        p=number(published.get("p_value",""),"p_value",meta_path,0,report)
        if p is not None and not 0<=p<=1:report.error("PROBABILITY_RANGE",test_id,meta_path)
    families=defaultdict(list)
    for rn,r in enumerate(meta_rows,2):
        p=number(r.get("p_value",""),"p_value",meta_path,rn,report);fdr=number(r.get("fdr",""),"fdr",meta_path,rn,report)
        if p is not None and fdr is not None:families[r.get("meta_family_id","")].append((rn,p,fdr))
    for family,items in families.items():
        for (rn,_,observed),expected in zip(items,bh([x[1] for x in items])):
            if not close_number(observed,expected):report.error("META_FDR_MISMATCH",family,meta_path,rn)

def validate_release(root,report,envelope_stage="final",allow_existing_envelope=False):
    required={"dataset_eligibility.tsv","published_results_audit.tsv","altitude_response_layers.tsv","cell_composition_sensitivity.tsv","altitude_pbmc_validation.tsv","whole_blood_hypotheses.tsv","release_status.json","artifact_index.tsv"}
    if envelope_stage=="final":required.add("_SUCCESS.json")
    for f in required:
        if not (root/f).is_file():report.error("MISSING_RELEASE_ARTIFACT",f,root/f)
    dirs=find_datasets(root,report); summaries={}
    for ds in EXPECTED_DATASETS:
        if ds not in dirs:report.error("MISSING_DATASET_DIRECTORY",ds,root)
        else:summaries[ds]=validate_dataset(dirs[ds],ds,report)
    elig=validate_eligibility(root/"dataset_eligibility.tsv",report)
    if set(dirs)!=set(elig):
        for ds in sorted(set(elig)-set(dirs)):report.error("ELIGIBILITY_WITHOUT_DATASET_DIRECTORY",ds,root/"dataset_eligibility.tsv")
        for ds in sorted(set(dirs)-set(elig)):report.error("DATASET_DIRECTORY_WITHOUT_ELIGIBILITY",ds,dirs[ds])
    for ds,d in dirs.items():
        if ds in elig:
            summary=summaries.get(ds,{});md=summary.get("metadata",{})
            try: es=int(elig[ds]["n_samples_included"]); ed=int(elig[ds]["n_donors_included"])
            except ValueError: continue
            if (es,ed)!=(md.get("n_samples"),md.get("n_donors")):report.error("ELIGIBILITY_COUNT_MISMATCH",ds,root/"dataset_eligibility.tsv")
            if elig[ds].get("dataset_manifest_sha256")!=summary.get("manifest_sha256") or elig[ds].get("dataset_payload_sha256")!=summary.get("dataset_payload_sha256"):report.error("ELIGIBILITY_DATASET_HASH_MISMATCH",ds,root/"dataset_eligibility.tsv")
            if elig[ds].get("analysis_role")!=summary.get("manifest",{}).get("analysis_role"):report.error("ELIGIBILITY_ROLE_MISMATCH",ds,root/"dataset_eligibility.tsv")
    def eligible_valid(ds):
        r=elig.get(ds,{});s=summaries.get(ds,{})
        return bool(s.get("valid")) and r.get("analysis_role")=="WHOLE_BLOOD_DE" and r.get("eligibility_status")=="ELIGIBLE" and str(r.get("de_authorized","")).lower() in TRUE and str(r.get("donor_independence_verified","")).lower() in TRUE and r.get("overlap_status")=="NO_OVERLAP"
    whole_blood_eligible=[ds for ds in ("GSE196728","GSE75665","GSE103940") if eligible_valid(ds)]
    status=load_json(root/"release_status.json",report) or {}
    declared=status.get("confirmation_dataset_ids",[])
    if not isinstance(declared,list):declared=[]
    valid_declared=(len(declared)==2 and len(set(declared))==2 and "GSE196728" in declared and len(set(declared)&CONFIRMATION_CANDIDATES)==1)
    if not valid_declared and declared:report.error("INVALID_CONFIRMATION_DATASET_SET",repr(declared),root/"release_status.json")
    eligible=declared if valid_declared and all(eligible_valid(ds) for ds in declared) else []
    phase=bool(eligible)
    report.metrics["phase_b_complete"]=phase; report.metrics["eligible_whole_blood_de_studies"]=whole_blood_eligible;report.metrics["confirmation_dataset_ids"]=eligible
    if bool(status.get("phase_b_complete"))!=phase:report.error("PHASE_B_STATUS_MISMATCH","",root/"release_status.json")
    if not phase:
        if status.get("phase_b_status")!="BLOCKED" or status.get("meta_status")!="NOT_RUN_INSUFFICIENT_STUDIES":report.error("PHASE_B_BLOCK_NOT_DECLARED","",root/"release_status.json")
        for name in ("altitude_whole_blood_meta.tsv","meta_study_effects.tsv","meta_method.json"):
            if (root/name).exists():report.error("META_FORBIDDEN_INSUFFICIENT_STUDIES",name,root/name)
    audit_rows,_=read_tsv(root/"published_results_audit.tsv",report,{"dataset_id","citation","contrast","reproduction_status","source_sha256","registry_version"})
    if {r.get("dataset_id") for r in audit_rows}!=set(dirs):report.error("AUDIT_DATASET_SET_MISMATCH","",root/"published_results_audit.tsv")
    for rn,r in enumerate(audit_rows,2):
        if r.get("dataset_id") not in dirs or not valid_hash(r.get("source_sha256","")) or any(missing(r.get(f)) for f in ("citation","contrast","reproduction_status","registry_version")):report.error("INVALID_AUDIT_ROW","",root/"published_results_audit.tsv",rn)
        elif r.get("source_sha256")!=sha(dirs[r["dataset_id"]]/"source_manifest.tsv"):report.error("AUDIT_SOURCE_HASH_MISMATCH","",root/"published_results_audit.tsv",rn)
    layer_rows,_=read_tsv(root/"altitude_response_layers.tsv",report,LAYER_REQ);layer_datasets={r.get("dataset_id") for r in layer_rows}
    if layer_datasets!=set(dirs):report.error("LAYER_DATASET_SET_MISMATCH","",root/"altitude_response_layers.tsv")
    for rn,r in enumerate(layer_rows,2):
        ds=r.get("dataset_id","")
        er=elig.get(ds,{})
        if r.get("analysis_role")!=summaries.get(ds,{}).get("manifest",{}).get("analysis_role") or r.get("eligibility_status")!=er.get("eligibility_status") or r.get("registry_version")!=er.get("registry_version"):report.error("LAYER_IDENTITY_FOREIGN_KEY_MISMATCH","",root/"altitude_response_layers.tsv",rn)
        if ds not in summaries or r.get("analysis_role") not in {"WHOLE_BLOOD_DE","PBMC_VALIDATION"} or r.get("analysis_state")!="VALIDATED":report.error("INVALID_LAYER_STATE","",root/"altitude_response_layers.tsv",rn)
        p=safe_child(root,r.get("output_path",""))
        if not p or not p.is_file() or sha(p)!=r.get("output_sha256"):report.error("LAYER_OUTPUT_HASH_MISMATCH","",root/"altitude_response_layers.tsv",rn)
        if r.get("output_path")!=f"{ds}/differential_expression.tsv":report.error("LAYER_OUTPUT_FOREIGN_KEY_MISMATCH","",root/"altitude_response_layers.tsv",rn)
        design_keys=summaries.get(ds,{}).get("design",{})
        if not any(k[1]==r.get("contrast_id") and k[2]==r.get("estimand") for k in design_keys):report.error("LAYER_DESIGN_FOREIGN_KEY_MISMATCH","",root/"altitude_response_layers.tsv",rn)
    for ds in set(dirs):
        expected_estimands={"pbmc_validation"} if ds=="GSE46480" else {"total_bulk","composition_conditional"}
        ds_rows=[r for r in layer_rows if r.get("dataset_id")==ds];observed={r.get("estimand") for r in ds_rows}
        if observed!=expected_estimands or len(ds_rows)!=len(expected_estimands):report.error("LAYER_ESTIMAND_SET_MISMATCH",f"{ds}: {sorted(observed)}",root/"altitude_response_layers.tsv")
    validate_composition(root,dirs,summaries,whole_blood_eligible,report)
    validate_pbmc(root,dirs,summaries,whole_blood_eligible,report)
    if phase:validate_meta(root,dirs,summaries,eligible,report)
    index_hash,payload=validate_artifact_index(root,report);report.metrics["payload_root_sha256"]=payload;report.metrics["envelope_stage"]=envelope_stage
    sealed_report=root/"validation_report.json";success_path=root/"_SUCCESS.json"
    if envelope_stage=="preseal":
        if not allow_existing_envelope and (sealed_report.exists() or success_path.exists()):report.error("PRESEAL_ENVELOPE_ALREADY_PRESENT","Remove sealed report and success marker before pre-seal validation",root)
        return
    if envelope_stage!="final":report.error("INVALID_ENVELOPE_STAGE",envelope_stage,root);return
    success=load_json(success_path,report) or {}
    sealer_path=Path(__file__).with_name("seal_module2.py")
    expected={"payload_root_sha256":payload,"artifact_index_sha256":index_hash,"validator_sha256":sha(Path(__file__)),"sealer_sha256":sha(sealer_path) if sealer_path.is_file() else "","schema_version":SCHEMA_VERSION,"schema_sha256":schema_sha256(),"validation_report_sha256":sha(sealed_report) if sealed_report.is_file() else "","validation_method":"CONTROLLED_PRESEAL_RERUN"}
    if not sealed_report.is_file():report.error("MISSING_SEALED_VALIDATION_REPORT","",sealed_report)
    else:
        sealed=load_json(sealed_report,report) or {}
        generated_at=sealed.get("generated_at_utc")
        if not valid_utc(generated_at):
            report.error("INVALID_SEALED_VALIDATION_REPORT","invalid generated_at_utc",sealed_report)
        else:
            replay=controlled_preseal_report(root,generated_at,allow_existing_envelope=True)
            if replay.error_count or sealed_report.read_bytes()!=report_bytes(replay.payload()):
                report.error("INVALID_SEALED_VALIDATION_REPORT","sealed bytes are not the complete current-validator pre-seal output",sealed_report)
    if set(success)!=set(expected) or any(success.get(k)!=v for k,v in expected.items()):report.error("SUCCESS_MARKER_MISMATCH","",success_path)
    else:report.metrics["success_marker_verified"]=True

def canonical_preseal_context(generated_at):
    return {"validation_stage":"preseal","report_outside_payload":True,"generated_at_utc":generated_at,"command":"release-preseal --release-dir <payload> --report <external>"}

def controlled_preseal_report(root,generated_at=None,allow_existing_envelope=False):
    replay=Report("release")
    if generated_at is not None:replay.generated_at_utc=generated_at
    replay.generation_context=canonical_preseal_context(replay.generated_at_utc)
    validate_release(root,replay,"preseal",allow_existing_envelope=allow_existing_envelope)
    return replay

def report_bytes(payload):
    return (json.dumps(payload,indent=2,sort_keys=True)+"\n").encode("utf-8")

def print_report(report):
    for d in report.diagnostics:
        loc=f" [{d.path}{':' + str(d.row) if d.row else ''}]" if d.path else ""; print(f"{d.severity} {d.code}: {d.message}{loc}")
    print(f"STATUS {'PASS' if not report.error_count else 'FAIL'}: {report.error_count} error(s), {report.warning_count} warning(s)")
def parser():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest="mode",required=True)
    x=s.add_parser("registry");x.add_argument("--project-root",type=Path,required=True);x.add_argument("--metadata",type=Path);x.add_argument("--report",type=Path)
    x=s.add_parser("dataset");x.add_argument("--dataset-dir",type=Path,required=True);x.add_argument("--dataset-id",required=True);x.add_argument("--report",type=Path)
    for name in ("release-preseal","release-final","release"):
        x=s.add_parser(name);x.add_argument("--release-dir",type=Path,required=True);x.add_argument("--report",type=Path,required=True)
    return p
def write_report(path,report):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+".tmp");tmp.write_bytes(report_bytes(report.payload()));tmp.replace(path)
def main(argv:Iterable[str]|None=None):
    a=parser().parse_args(argv);release_mode=a.mode in {"release-preseal","release-final","release"};r=Report("release" if release_mode else a.mode)
    if a.mode=="registry":validate_registry(a.project_root.resolve(),a.metadata.resolve() if a.metadata else None,r)
    elif a.mode=="dataset":validate_dataset(a.dataset_dir.resolve(),a.dataset_id,r)
    else:
        root=a.release_dir.resolve();out=a.report.resolve()
        try:
            out.relative_to(root);r.error("REPORT_INSIDE_PAYLOAD","--report must be outside release payload",out);print_report(r);return 1
        except ValueError:pass
        stage="preseal" if a.mode=="release-preseal" else "final"
        r.generation_context=canonical_preseal_context(r.generated_at_utc) if stage=="preseal" else {"validation_stage":"final","report_outside_payload":True,"generated_at_utc":r.generated_at_utc,"command":"release-final --release-dir <payload> --report <external>"}
        validate_release(root,r,stage)
    if getattr(a,"report",None):write_report(a.report.resolve(),r)
    print_report(r);return 0 if not r.error_count else 1
if __name__=="__main__":sys.exit(main())
