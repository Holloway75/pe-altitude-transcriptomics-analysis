#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, importlib.util, json, math, shutil, sys, tempfile, unittest
from unittest import mock
from pathlib import Path

# Moved from .claude/skills/module2-geo-expression/tests/ into the project on
# 2026-09-03; the validator now sits next to this test instead of one level up.
SCRIPT=Path(__file__).resolve().with_name("validate_module2.py")
SPEC=importlib.util.spec_from_file_location("module2_validator",SCRIPT);assert SPEC and SPEC.loader
V=importlib.util.module_from_spec(SPEC);sys.modules[SPEC.name]=V;SPEC.loader.exec_module(V)
PUBLISH_SCRIPT=SCRIPT.with_name("seal_module2.py")
PUBLISH_SPEC=importlib.util.spec_from_file_location("module2_sealer",PUBLISH_SCRIPT);assert PUBLISH_SPEC and PUBLISH_SPEC.loader
P=importlib.util.module_from_spec(PUBLISH_SPEC);sys.modules[PUBLISH_SPEC.name]=P;PUBLISH_SPEC.loader.exec_module(P)
H="1"*64

def tsv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def rows(path):
    with path.open(encoding="utf-8",newline="") as h:return list(csv.DictReader(h,delimiter="\t"))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def file_hash_map(d):return {p.relative_to(d).as_posix():digest(p) for p in d.rglob("*") if p.is_file() and p.name!="dataset_manifest.json"}
def write_dataset_manifest(d,dataset,role):
    hashes=file_hash_map(d);payload=hashlib.sha256("".join(f"{k}\t{hashes[k]}\n" for k in sorted(hashes)).encode()).hexdigest()
    m={"schema_version":V.SCHEMA_VERSION,"dataset_id":dataset,"analysis_role":role,"registry_version":"r1","canonical_directory":dataset,"artifact_hashes":hashes,"dataset_payload_sha256":payload};(d/"dataset_manifest.json").write_text(json.dumps(m));return m

def metadata(dataset,donors=2):
    if dataset=="GSE196728":cells=[(p,c) for p in ("sea_level","3800m","5100m") for c in V.GSE196728_CLOCKS]
    elif dataset=="GSE46480":cells=[("baseline","NA"),("day3","NA")]
    elif dataset=="GSE103940":cells=[("plain","NA"),("highaltitude","NA")]
    else:cells=[("plain","NA"),("5300m","NA")]
    out=[]
    for di in range(1,donors+1):
        for i,(level,clock) in enumerate(cells):
            out.append({"dataset_id":dataset,"sample_id":f"{dataset}_S{di}_{i}","expression_column":f"E{di}_{i}","donor_id":f"D{di}","include":"true","exclusion_reason":"NA","tissue":"PBMC" if dataset=="GSE46480" else "Whole Blood","assay":"microarray" if dataset=="GSE46480" else "RNA-seq","altitude_label":level,"pair_id":f"D{di}","source_url":"https://ncbi.example/source","planned_imbalance_id":"NA","acetazolamide_raw":"no" if dataset=="GSE46480" else "NA","time_of_day":clock})
    return out

def design_files(d,md,estimands,engine="fixed_effect"):
    samples=[r["sample_id"] for r in md];donor_order=[r["donor_id"] for r in md];donors=sorted(set(donor_order));ref="baseline" if md[0]["dataset_id"]=="GSE46480" else "sea_level" if md[0]["dataset_id"]=="GSE196728" else "plain"
    cols=["intercept"]+([] if engine!="fixed_effect" else [f"donor_{x}" for x in donors[1:]])+["condition"]
    matrix=[]
    for r in md:matrix.append({"sample_id":r["sample_id"],**dict(zip(cols,[1]+([] if engine!="fixed_effect" else [int(r["donor_id"]==x) for x in donors[1:]])+[int(r["altitude_label"]!=ref)]))})
    tsv(d/"design_matrix.tsv",matrix);tsv(d/"contrast.tsv",[{"contrast":"target",**{c:int(c=="condition") for c in cols}}]);tsv(d/"block_vector.tsv",[{"sample_id":s,"donor_id":x} for s,x in zip(samples,donor_order)])
    rank=V.matrix_rank([[r[c] for c in cols] for r in matrix]);sample_hash=hashlib.sha256(("\n".join(samples)+"\n").encode()).hexdigest();donor_hash=hashlib.sha256(("\n".join(donor_order)+"\n").encode()).hexdigest();out=[]
    for est in estimands:
        out.append({"dataset_id":md[0]["dataset_id"],"model_id":f"m_{est}","estimand":est,"formula":"~ donor + condition" if engine=="fixed_effect" else "~ condition","contrast_id":"target","contrast_definition":"condition","scientific_comparison_id":"altitude_primary","reference_level":ref,"n_samples":len(md),"n_donors":len(donors),"design_rank":rank,"design_columns":len(cols),"residual_df":len(md)-rank,"correlation_engine":engine,"correlation_structure":"donor fixed" if engine=="fixed_effect" else "donor block","block_vector_path":"NA" if engine=="fixed_effect" else "block_vector.tsv","block_vector_sha256":"NA" if engine=="fixed_effect" else digest(d/"block_vector.tsv"),"donor_key_sha256":donor_hash,"random_effect_formula":"(1|donor_id)" if engine=="dream" else "NA","primary":"true","interpretation_boundary":"location order acclimatization" if md[0]["dataset_id"]=="GSE196728" else "paired study","registry_version":"r1","design_matrix_path":"design_matrix.tsv","design_matrix_sha256":digest(d/"design_matrix.tsv"),"sample_order_sha256":sample_hash,"contrast_vector_path":"contrast.tsv","contrast_vector_sha256":digest(d/"contrast.tsv")})
    tsv(d/"design.tsv",out);return out

def build_dataset(root,dataset="GSE46480",donors=2,engine="fixed_effect"):
    d=root/dataset;md=metadata(dataset,donors);tsv(d/"sample_metadata.tsv",md);estimands=["total_bulk","composition_conditional"] if dataset in {"GSE196728","GSE75665","GSE103940"} else ["pbmc_validation"];design_files(d,md,estimands,engine)
    features=[{"dataset_id":dataset,"feature_id":f"F{i}","feature_namespace":"Ensembl","feature_namespace_version":"v1"} for i in range(4)];tsv(d/"feature_manifest.tsv",features);feature_ids=[x["feature_id"] for x in features]
    de=[];native=[];universes=[];families=[];ps=[.01,.04,.2];qs=[.03,.06,.2];effects=[.8,-.3,.1] if dataset=="GSE103940" else [1,-.5,0]
    for est in estimands:
        fid=f"fam_{est}";family_rows=[]
        for i,feature in enumerate(feature_ids):
            tested=i<3;test=f"T{i}";status="TESTED" if tested else "FILTERED";family_rows.append((feature,test,status))
        universe_hash=hashlib.sha256("".join(f"{f}\t{t}\t{s}\n" for f,t,s in family_rows).encode()).hexdigest();families.append({"bh_family_id":fid,"model_id":f"m_{est}","contrast_id":"target","estimand":est,"test_universe_sha256":universe_hash})
        for i,(feature,test,status) in enumerate(family_rows):
            tested=status=="TESTED";e=effects[i] if tested else "NA";p=ps[i] if tested else "NA";q=qs[i] if tested else "NA";stat=e/.2 if tested else "NA"
            de.append({"dataset_id":dataset,"model_id":f"m_{est}","contrast_id":"target","estimand":est,"test_id":test,"feature_id":feature,"feature_namespace":"Ensembl","feature_namespace_version":"v1","gene_id":f"G{i}","gene_symbol":f"SYM{i}","mapping_status":"one_to_one","log2fc":e,"se":.2 if tested else "NA","statistic":stat,"statistic_type":"moderated_t","df":10 if tested else "NA","p_value":p,"fdr":q,"direction":"up" if tested and e>0 else "down" if tested and e<0 else "zero","deg":"true" if tested and q<=.05 and abs(e)>=.5 else "false","base_expression":5 if tested else 0,"n_samples":len(md),"n_donors":donors,"annotation_version":"a1","method_version":"m1","test_status":status,"not_tested_reason":"NA" if tested else "PREDECLARED_LOW_EXPRESSION","bh_family_id":fid,"test_universe_sha256":universe_hash})
            if tested:native.append({"dataset_id":dataset,"model_id":f"m_{est}","contrast_id":"target","estimand":est,"test_id":test,"feature_id":feature,"log2fc":e,"statistic":stat,"df":10,"p_value":p})
        tids=[f"T{i}" for i in range(3)];universes.append({"dataset_id":dataset,"model_id":f"m_{est}","contrast_id":"target","estimand":est,"bh_family_id":fid,"input_feature_count":4,"tested_count":3,"input_feature_sha256":digest(d/"feature_manifest.tsv"),"tested_feature_sha256":hashlib.sha256(("\n".join(tids)+"\n").encode()).hexdigest(),"family_universe_sha256":universe_hash})
    tsv(d/"differential_expression.tsv",de);tsv(d/"native_results.tsv",native);tsv(d/"test_universe.tsv",universes);(d/"raw.dat").write_bytes(b"x");(d/"minimum_gate.md").write_text("Synthetic precision/power gate.\n")
    lock={"R":{"Version":"4.5.0"},"Bioconductor":{"Version":"3.21"},"Repositories":[{"Name":"CRAN","URL":"snapshot"},{"Name":"BioC","URL":"bioc-snapshot"}],"Packages":{"limma":{"Package":"limma","Version":"3.0","Source":"Bioconductor","Imports":["BiocGenerics"],"LinkingTo":[]},"BiocGenerics":{"Package":"BiocGenerics","Version":"0.1","Source":"Bioconductor","Depends":["methods"]}}};(d/"renv.lock").write_text(json.dumps(lock))
    source={"dataset_id":dataset,"source_kind":"GEO","url":"https://ncbi.example/source","retrieved_at_utc":"2026-07-23T00:00:00Z","local_path":"raw.dat","bytes":1,"sha256":digest(d/"raw.dat"),"content_role":"normalized_expression" if dataset=="GSE46480" else "counts","archive_check":"NOT_APPLICABLE","access_status":"PUBLIC","record_status_at_utc":"2026-07-23T00:00:00Z","citation":"PMID:1","terms_url":"https://ncbi.example/terms","license_status":"NOT_STATED","use_restrictions":"repository terms","privacy_handling":"no re-identification","local_verification_status":"VERIFIED","remote_immutable_id":dataset,"verification_reason":"local hash verified","audit_exception_status":"NOT_REQUIRED","audit_exception_id":"NA","audit_exception_approved_by":"NA"};tsv(d/"source_manifest.tsv",[source])
    expected_samples,expected_donors=V.EXPECTED_FLOW[dataset];quant={"source_type":"MICROARRAY_LOG2","quantifier":"array","tx2gene_sha256":H,"counts_from_abundance":"not_applicable","offset_decision":"not_applicable","constructor":"ExpressionSet"} if dataset=="GSE46480" else {"source_type":"INTEGER_COUNTS","quantifier":"gene_counts","tx2gene_sha256":H,"counts_from_abundance":"not_applicable","offset_decision":"not_applicable","constructor":"DGEList_voom"}
    method={"schema_version":V.SCHEMA_VERSION,"dataset_id":dataset,"registry_version":"r1","method_version":"m1","engine":"limma" if dataset=="GSE46480" else "limma_voom","correlation_engine":engine,"r_version":"4.5.0","bioconductor_version":"3.21","repositories":{"CRAN":"snapshot","BioC":"bioc-snapshot"},"environment_lock":{"type":"renv","path":"renv.lock","sha256":digest(d/"renv.lock")},"packages":{"limma":"3.0"},"session_info":"R session","bioc_valid":True,"normalization":{"method":"frozen"},"filter":{"method":"filterByExpr","outcome_independent":True},"annotation":{"version":"a1","sha256":H},"composition":{"method":"frozen"},"quantification":quant,"geoquery":{"package_version":"2.80.0","calls":[{"function":"getGEO","GSEMatrix":True}],"cache_manifest_sha256":H,"platforms":["GPLTEST"],"alignment_verified":True,"supplement_classification":"content_inspected"},"contrasts":["target"],"seeds":{"composition":1},"deg_rule":{"fdr_lte":.05,"abs_log2fc_gte":.5},"test_families":families,"native_results":[{"statistic_type":"moderated_t","path":"native_results.tsv","sha256":digest(d/"native_results.tsv")}],"input_hashes":{"raw.dat":digest(d/"raw.dat"),"feature_manifest.tsv":digest(d/"feature_manifest.tsv")},"output_hashes":{"differential_expression.tsv":digest(d/"differential_expression.tsv"),"test_universe.tsv":digest(d/"test_universe.tsv")},"sample_flow":{"expected_samples":expected_samples,"expected_donors":expected_donors,"minimum_donors":2,"maximum_attrition":expected_samples-len(md),"exclusion_reasons_complete":True,"minimum_gate_rationale_path":"minimum_gate.md","minimum_gate_rationale_sha256":digest(d/"minimum_gate.md")}};(d/"method.json").write_text(json.dumps(method))
    design_rows=rows(d/"design.tsv");qc_values={"included_samples":len(md),"included_donors":donors,"excluded_samples":0,"input_feature_count":4,"tested_count_total":3*len(estimands),"design_rank_min":min(int(x["design_rank"]) for x in design_rows),"residual_df_min":min(int(x["residual_df"]) for x in design_rows),"test_universe_sha256":digest(d/"test_universe.tsv")};qc=[{"dataset_id":dataset,"stage":"gate","metric":k,"value":v,"status":"PASS","threshold":"frozen","details":"synthetic"} for k,v in qc_values.items()];qc.extend([{"dataset_id":dataset,"stage":"gate","metric":k,"value":1,"status":"PASS","threshold":"frozen","details":"synthetic"} for k in ("exclusion_reasons_complete","design_rank","residual_df","normalization","mapping_coverage")]);tsv(d/"qc_summary.tsv",qc)
    role="WHOLE_BLOOD_DE" if dataset in {"GSE196728","GSE75665","GSE103940"} else "PBMC_VALIDATION";write_dataset_manifest(d,dataset,role);return d

def add_second_total_bulk_family(d):
    design=rows(d/"design.tsv");base=next(r for r in design if r["estimand"]=="total_bulk");second=dict(base);second.update({"model_id":"m_total_bulk_secondary","contrast_id":"target_secondary","scientific_comparison_id":"altitude_secondary","primary":"false"});design.append(second);tsv(d/"design.tsv",design)
    de=rows(d/"differential_expression.tsv");copies=[]
    for r in [x for x in de if x["estimand"]=="total_bulk"]:
        q=dict(r);q.update({"model_id":"m_total_bulk_secondary","contrast_id":"target_secondary","bh_family_id":"fam_total_bulk_secondary"})
        if q["test_status"]=="TESTED":
            q["log2fc"]=str(-float(q["log2fc"]));q["statistic"]=str(-float(q["statistic"]));q["direction"]="up" if float(q["log2fc"])>0 else "down" if float(q["log2fc"])<0 else "zero"
        copies.append(q)
    de.extend(copies);tsv(d/"differential_expression.tsv",de)
    native=rows(d/"native_results.tsv")
    for r in [x for x in native if x["model_id"]=="m_total_bulk"]:
        q=dict(r);q.update({"model_id":"m_total_bulk_secondary","contrast_id":"target_secondary","log2fc":str(-float(q["log2fc"])),"statistic":str(-float(q["statistic"]))});native.append(q)
    tsv(d/"native_results.tsv",native)
    universe=rows(d/"test_universe.tsv");u=dict(next(r for r in universe if r["bh_family_id"]=="fam_total_bulk"));u.update({"model_id":"m_total_bulk_secondary","contrast_id":"target_secondary","bh_family_id":"fam_total_bulk_secondary"});universe.append(u);tsv(d/"test_universe.tsv",universe)
    method=json.loads((d/"method.json").read_text());family=dict(next(r for r in method["test_families"] if r["bh_family_id"]=="fam_total_bulk"));family.update({"bh_family_id":"fam_total_bulk_secondary","model_id":"m_total_bulk_secondary","contrast_id":"target_secondary"});method["test_families"].append(family);method["contrasts"].append("target_secondary");method["native_results"][0]["sha256"]=digest(d/"native_results.tsv");method["output_hashes"].update({"differential_expression.tsv":digest(d/"differential_expression.tsv"),"test_universe.tsv":digest(d/"test_universe.tsv")});(d/"method.json").write_text(json.dumps(method))
    qc=rows(d/"qc_summary.tsv")
    for r in qc:
        if r["metric"]=="tested_count_total":r["value"]=str(int(r["value"])+3)
        if r["metric"]=="test_universe_sha256":r["value"]=digest(d/"test_universe.tsv")
    tsv(d/"qc_summary.tsv",qc);write_dataset_manifest(d,d.name,"WHOLE_BLOOD_DE")

def elig_row(d,dataset,status="ELIGIBLE"):
    m=json.loads((d/"dataset_manifest.json").read_text());md=rows(d/"sample_metadata.tsv");role=m["analysis_role"];expected=V.EXPECTED_FLOW[dataset];donors=len({r["donor_id"] for r in md if r["include"]=="true"})
    return {"dataset_id":dataset,"evidence_layer":"layer","eligibility_status":status,"n_samples_included":len(md),"n_donors_included":donors,"pairing_complete":"true","donor_independence_verified":"true","tissue_match":"true","primary_contrast_estimable":"true","design_full_rank":"true","residual_df_adequate":"true","source_data_available":"true","overlap_status":"NO_OVERLAP","decision_reason":"frozen","registry_version":"r1","analysis_role":role,"expected_samples":expected[0],"expected_donors":expected[1],"minimum_donors":2,"maximum_attrition":expected[0]-len(md),"exclusion_count":0,"exclusion_reasons_complete":"true","de_authorized":"true","no_de_reason":"NA","dataset_manifest_sha256":digest(d/"dataset_manifest.json"),"dataset_payload_sha256":m["dataset_payload_sha256"]}

def build_index(root):
    for p in (root/"artifact_index.tsv",root/"_SUCCESS.json",root/"validation_report.json"):
        if p.exists():p.unlink()
    payload=[p for p in root.rglob("*") if p.is_file() and p.name not in {"artifact_index.tsv","_SUCCESS.json","validation_report.json"}];index=[{"path":p.relative_to(root).as_posix(),"bytes":p.stat().st_size,"sha256":digest(p)} for p in payload];tsv(root/"artifact_index.tsv",index)
    return hashlib.sha256(V.canonical_index(index).encode()).hexdigest()

def marker(root,payload_hash):
    success={"payload_root_sha256":payload_hash,"artifact_index_sha256":digest(root/"artifact_index.tsv"),"validator_sha256":digest(SCRIPT),"schema_version":V.SCHEMA_VERSION,"schema_sha256":V.schema_sha256(),"validation_report_sha256":digest(root/"validation_report.json")};(root/"_SUCCESS.json").write_text(json.dumps(success))

def publisher_seal(root):
    payload_hash=build_index(root);outside=root.parent/f"{root.name}.preseal.json";r=V.Report("release");r.generation_context={"validation_stage":"preseal","report_outside_payload":True,"generated_at_utc":r.generated_at_utc,"command":"release-preseal --release-dir <payload> --report <external>"};V.validate_release(root,r,"preseal")
    if r.error_count:raise AssertionError(r.payload())
    V.write_report(outside,r)
    if P.main(["--release-dir",str(root),"--preseal-report",str(outside)])!=0:raise AssertionError("publisher rejected valid pre-seal report")
    return outside

def forge_seal(root):
    payload_hash=build_index(root);sealed={"schema_version":V.SCHEMA_VERSION,"schema_sha256":V.schema_sha256(),"validator_sha256":digest(SCRIPT),"mode":"release","status":"PASS","schema_pass":True,"scientific_gate_pass":True,"payload_root_sha256":payload_hash,"generation_context":{"validation_stage":"preseal","report_outside_payload":True,"generated_at_utc":"2026-07-23T00:00:00Z","command":"forged test envelope"}};(root/"validation_report.json").write_text(json.dumps(sealed));marker(root,payload_hash)

def build_release(root,phase_complete=False,sealed=True):
    d196=build_dataset(root,"GSE196728",2);d756=build_dataset(root,"GSE75665",2);d103=build_dataset(root,"GSE103940",2);d464=build_dataset(root,"GSE46480",2);dirs=[d196,d756,d103,d464]
    tsv(root/"dataset_eligibility.tsv",[elig_row(d,d.name) for d in dirs]);tsv(root/"published_results_audit.tsv",[{"dataset_id":d.name,"citation":"PMID","contrast":"target","reproduction_status":"FROZEN","source_sha256":digest(d/"source_manifest.tsv"),"registry_version":"r1"} for d in dirs])
    layers=[]
    for d in dirs:
        ds=d.name
        role="PBMC_VALIDATION" if ds=="GSE46480" else "WHOLE_BLOOD_DE";ests=["pbmc_validation"] if ds=="GSE46480" else ["total_bulk","composition_conditional"]
        for est in ests:layers.append({"dataset_id":ds,"analysis_role":role,"analysis_state":"VALIDATED","contrast_id":"target","estimand":est,"output_path":f"{ds}/differential_expression.tsv","output_sha256":digest(d/"differential_expression.tsv"),"eligibility_status":"ELIGIBLE","registry_version":"r1"})
    tsv(root/"altitude_response_layers.tsv",layers)
    whole=[d196,d756,d103];comp=[]
    for d in whole:
        de=rows(d/"differential_expression.tsv");by={(r["estimand"],r["test_id"]):r for r in de if r["test_status"]=="TESTED"}
        for i in range(3):
            total=by[("total_bulk",f"T{i}")];conditional=by[("composition_conditional",f"T{i}")];delta,attenuation,relation,sensitivity=V.composition_values(float(total["log2fc"]),float(conditional["log2fc"]))
            comp.append({"dataset_id":d.name,"total_model_id":"m_total_bulk","conditional_model_id":"m_composition_conditional","contrast_id":"target","test_id":f"T{i}","total_bulk_log2fc":total["log2fc"],"total_bulk_se":total["se"],"total_bulk_statistic":total["statistic"],"total_bulk_fdr":total["fdr"],"conditional_log2fc":conditional["log2fc"],"conditional_se":conditional["se"],"conditional_statistic":conditional["statistic"],"conditional_fdr":conditional["fdr"],"delta_log2fc":delta,"attenuation_fraction":"NA" if attenuation is None else attenuation,"direction_relation":relation,"sensitivity_class":sensitivity,"source_de_sha256":digest(d/"differential_expression.tsv")})
    tsv(root/"cell_composition_sensitivity.tsv",comp)
    hypotheses=[{"target_id":"PROGRAM1","registry_version":"r1","frozen_at_utc":"2026-07-23T00:00:00Z","source_dataset_id":"GSE196728","source_model_id":"m_total_bulk","source_contrast_id":"target","source_estimand":"total_bulk","source_test_id":"T0","source_de_sha256":digest(d196/"differential_expression.tsv"),"hypothesis_family_id":"wb_programs","pbmc_model_id":"m_pbmc_validation","pbmc_contrast_id":"target","pbmc_estimand":"pbmc_validation","pbmc_test_id":"T0","multiplicity_family_id":"fam_pbmc_validation"}];tsv(root/"whole_blood_hypotheses.tsv",hypotheses)
    pbmc=next(r for r in rows(d464/"differential_expression.tsv") if r["estimand"]=="pbmc_validation" and r["test_id"]=="T0")
    tsv(root/"altitude_pbmc_validation.tsv",[{"dataset_id":"GSE46480","model_id":"m_pbmc_validation","contrast_id":"target","estimand":"pbmc_validation","test_id":"T0","target_id":"PROGRAM1","target_source_sha256":digest(root/"whole_blood_hypotheses.tsv"),"hypothesis_status":"PRESPECIFIED_TARGET","statistic":pbmc["statistic"],"p_value":pbmc["p_value"],"fdr":pbmc["fdr"],"multiplicity_family_id":"fam_pbmc_validation","registry_version":"r1"}])
    phase=phase_complete;(root/"release_status.json").write_text(json.dumps({"phase_b_complete":phase,"phase_b_status":"COMPLETE" if phase else "BLOCKED","meta_status":"COMPLETE" if phase else "NOT_RUN_INSUFFICIENT_STUDIES","confirmation_dataset_ids":["GSE196728","GSE103940"] if phase else []}))
    if phase_complete:
        eligible=[d196,root/"GSE103940"];contributors=[];effects=[]
        for d in eligible:
            manifest=json.loads((d/"dataset_manifest.json").read_text());method=json.loads((d/"method.json").read_text());family=next(x for x in method["test_families"] if x["bh_family_id"]=="fam_total_bulk");contributors.append({"dataset_id":d.name,"source_family_id":"fam_total_bulk","model_id":"m_total_bulk","contrast_id":"target","estimand":"total_bulk","scientific_comparison_id":"altitude_primary","test_universe_sha256":family["test_universe_sha256"],"source_de_sha256":digest(d/"differential_expression.tsv"),"dataset_payload_sha256":manifest["dataset_payload_sha256"]})
            for source in rows(d/"differential_expression.tsv"):
                if source["estimand"]=="total_bulk" and source["test_status"]=="TESTED":effects.append({"dataset_id":d.name,"source_family_id":"fam_total_bulk","model_id":"m_total_bulk","contrast_id":"target","scientific_comparison_id":"altitude_primary","test_id":source["test_id"],"estimand":"total_bulk","effect":source["log2fc"],"se":source["se"],"source_de_sha256":digest(d/"differential_expression.tsv")})
        tsv(root/"meta_study_effects.tsv",effects);native=[]
        for i in range(3):
            test=f"T{i}";source=[r for r in effects if r["test_id"]==test];stats=V.meta_statistics(source);native.append({"meta_family_id":"altitude_primary_total_bulk","scientific_comparison_id":"altitude_primary","estimand":"total_bulk","test_id":test,"k_studies":2,**stats,"fdr":0,"estimand_comparability_status":"COMPARABLE","source_effects_sha256":V.source_effects_sha256(source)})
        adjusted=V.bh([r["p_value"] for r in native])
        for r,q in zip(native,adjusted):r["fdr"]=q
        tsv(root/"meta_native_results.tsv",native);tsv(root/"altitude_whole_blood_meta.tsv",native)
        (root/"meta_method.json").write_text(json.dumps({"k":2,"tau2_estimator":"REML","interval_method":"HKSJ","interpretation":"DESCRIPTIVE","confirmatory_deg":False,"sensitivity_methods":["modified_HKSJ","fixed_effect"],"meta_family_id":"altitude_primary_total_bulk","scientific_comparison_id":"altitude_primary","estimand":"total_bulk","contributors":contributors,"native_results":{"path":"meta_native_results.tsv","sha256":digest(root/"meta_native_results.tsv"),"package":"metafor","package_version":"4.8-0","statistic_contract":"REML_HKSJ"}}))
    if sealed:publisher_seal(root)
    else:build_index(root)
    return root

class Tests(unittest.TestCase):
    def dataset_report(self,d,ds):r=V.Report("dataset");V.validate_dataset(d,ds,r);return r
    def release_report(self,root):r=V.Report("release");V.validate_release(root,r);return r
    def test_valid_dataset_and_blocked_release(self):
        with tempfile.TemporaryDirectory() as x:
            self.assertEqual(self.dataset_report(build_dataset(Path(x)),"GSE46480").error_count,0);r=self.release_report(build_release(Path(x)/"r"));self.assertEqual(r.error_count,0,r.payload());self.assertFalse(r.metrics["phase_b_complete"])
    def test_frozen_confirmation_pair_unlocks_phase_b(self):
        with tempfile.TemporaryDirectory() as x:
            r=self.release_report(build_release(Path(x),True));self.assertEqual(r.error_count,0,r.payload());self.assertTrue(r.metrics["phase_b_complete"])
    def test_registry_infers_dataset_and_rejects_bad_grid(self):
        with tempfile.TemporaryDirectory() as x:
            root=Path(x);tsv(root/"config/expression_datasets.tsv",[{"dataset_id":ds,"tissue":v["tissue"],"design_role":"x","analysis_role":v["role"],"de_authorized":"true"} for ds,v in V.EXPECTED_DATASETS.items()]);bad=metadata("GSE75665",10)
            for r in bad:r["altitude_label"]="badA" if r["altitude_label"]=="plain" else "badB"
            tsv(root/"bad.tsv",bad);r=V.Report("registry");V.validate_registry(root,root/"bad.tsv",r);self.assertIn("GSE75665_PAIR_GRID",{d.code for d in r.diagnostics})
    def test_whole_clock_level_missing_requires_imbalance(self):
        with tempfile.TemporaryDirectory() as x:
            p=Path(x)/"m.tsv";md=metadata("GSE196728",2);md=[r for r in md if not (r["donor_id"]=="D2" and r["time_of_day"]=="22")];tsv(p,md);r=V.Report("metadata");V.validate_metadata(p,r,"GSE196728");self.assertIn("UNDECLARED_IMBALANCE",{d.code for d in r.diagnostics})
    def test_gse196728_uses_primary_metadata_local_clocks(self):
        self.assertEqual(V.GSE196728_CLOCKS,("02","06","10","14","18","22"))
        with tempfile.TemporaryDirectory() as x:
            p=Path(x)/"m.tsv";md=metadata("GSE196728",2)
            shifted=dict(zip(V.GSE196728_CLOCKS,("00","04","08","12","16","20")))
            for r in md:r["time_of_day"]=shifted[r["time_of_day"]]
            tsv(p,md);report=V.Report("metadata");V.validate_metadata(p,report,"GSE196728")
            self.assertIn("GSE196728_CLOCK_SET",{d.code for d in report.diagnostics})
    def test_duplicate_correlation_branch_passes_and_bad_block_fails(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x),"GSE196728",2,"duplicateCorrelation");self.assertEqual(self.dataset_report(d,"GSE196728").error_count,0);br=rows(d/"block_vector.tsv");br[0]["donor_id"]="WRONG";tsv(d/"block_vector.tsv",br);dr=rows(d/"design.tsv");
            for r in dr:r["block_vector_sha256"]=digest(d/"block_vector.tsv")
            tsv(d/"design.tsv",dr);self.assertIn("BLOCK_VECTOR_MISMATCH",{q.code for q in self.dataset_report(d,"GSE196728").diagnostics})
    def test_mixed_fixed_and_block_evidence_fails(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x));dr=rows(d/"design.tsv");dr[0]["block_vector_path"]="block_vector.tsv";dr[0]["block_vector_sha256"]=digest(d/"block_vector.tsv");tsv(d/"design.tsv",dr);self.assertIn("CONFLICTING_CORRELATION_EVIDENCE",{q.code for q in self.dataset_report(d,"GSE46480").diagnostics})
    def test_dream_branch_passes_and_invalid_formula_fails(self):
        with tempfile.TemporaryDirectory() as x:
            root=Path(x);d=build_dataset(root,"GSE196728",2,"dream");self.assertEqual(self.dataset_report(d,"GSE196728").error_count,0);dr=rows(d/"design.tsv");dr[0]["random_effect_formula"]="(1|wrong_id)";tsv(d/"design.tsv",dr);self.assertIn("INVALID_RANDOM_EFFECT_FORMULA",{q.code for q in self.dataset_report(d,"GSE196728").diagnostics});d2=build_dataset(root/"mismatch","GSE196728",2,"dream");m=json.loads((d2/"method.json").read_text());m["correlation_engine"]="fixed_effect";(d2/"method.json").write_text(json.dumps(m));self.assertIn("METHOD_CORRELATION_ENGINE_MISMATCH",{q.code for q in self.dataset_report(d2,"GSE196728").diagnostics})
    def test_truncated_feature_universe_and_duplicate_feature_fail(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x));de=rows(d/"differential_expression.tsv");de=de[:-1];de[1]["feature_id"]=de[0]["feature_id"];tsv(d/"differential_expression.tsv",de);codes={q.code for q in self.dataset_report(d,"GSE46480").diagnostics};self.assertIn("INCOMPLETE_FEATURE_UNIVERSE",codes);self.assertIn("DUPLICATE_FEATURE_IN_FAMILY",codes)
    def test_native_p_and_dataset_id_mismatch_fail(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x));de=rows(d/"differential_expression.tsv");de[0]["p_value"]=".5";de[0]["dataset_id"]="GSE75665";tsv(d/"differential_expression.tsv",de);codes={q.code for q in self.dataset_report(d,"GSE46480").diagnostics};self.assertIn("NATIVE_RESULT_VALUE_MISMATCH",codes);self.assertIn("DE_DATASET_ID_MISMATCH",codes)
    def test_method_output_hash_qc_and_version_foreign_keys_fail(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x));m=json.loads((d/"method.json").read_text());m["output_hashes"]["differential_expression.tsv"]=H;(d/"method.json").write_text(json.dumps(m));de=rows(d/"differential_expression.tsv");de[0]["annotation_version"]="wrong";tsv(d/"differential_expression.tsv",de);q=rows(d/"qc_summary.tsv");q[0]["value"]="999";tsv(d/"qc_summary.tsv",q);codes={z.code for z in self.dataset_report(d,"GSE46480").diagnostics};self.assertTrue({"METHOD_FILE_HASH_MISMATCH","ANNOTATION_VERSION_MISMATCH","QC_VALUE_MISMATCH"}.issubset(codes))
    def test_quantification_allow_matrix_rejects_continuous_count(self):
        with tempfile.TemporaryDirectory() as x:
            root=Path(x);d=build_dataset(root);m=json.loads((d/"method.json").read_text());m["engine"]="DESeq2";m["quantification"].update({"source_type":"NORMALIZED_CONTINUOUS","constructor":"DESeqDataSetFromMatrix"});(d/"method.json").write_text(json.dumps(m));self.assertIn("INVALID_QUANTIFICATION_COMBINATION",{q.code for q in self.dataset_report(d,"GSE46480").diagnostics});d2=build_dataset(root/"bad_quantifier");m=json.loads((d2/"method.json").read_text());m["quantification"]["quantifier"]="salmon";(d2/"method.json").write_text(json.dumps(m));self.assertIn("INVALID_QUANTIFICATION_COMBINATION",{q.code for q in self.dataset_report(d2,"GSE46480").diagnostics})
    def test_environment_lock_cross_checks_versions_dependencies(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x));lock=json.loads((d/"renv.lock").read_text());lock["R"]["Version"]="4.4";lock["Packages"]["limma"]["Imports"]=["missingPkg"];(d/"renv.lock").write_text(json.dumps(lock));m=json.loads((d/"method.json").read_text());m["environment_lock"]["sha256"]=digest(d/"renv.lock");(d/"method.json").write_text(json.dumps(m));codes={q.code for q in self.dataset_report(d,"GSE46480").diagnostics};self.assertIn("LOCK_R_VERSION_MISMATCH",codes);self.assertIn("LOCK_DEPENDENCY_MISSING",codes)
    def test_unknown_access_and_unverified_remote_require_approval(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x));s=rows(d/"source_manifest.tsv");s[0].update({"access_status":"UNKNOWN","local_verification_status":"REMOTE_ONLY","use_restrictions":"NA","audit_exception_status":"NOT_APPROVED"});tsv(d/"source_manifest.tsv",s);codes={q.code for q in self.dataset_report(d,"GSE46480").diagnostics};self.assertIn("SOURCE_AUDIT_APPROVAL_REQUIRED",codes);self.assertIn("INCOMPLETE_GOVERNANCE",codes)
    def test_duplicate_dataset_directory_and_orphan_artifact_fail(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x));copy=root/"nested"/"GSE75665";copy.parent.mkdir();shutil.copytree(root/"GSE75665",copy);(root/"orphan").mkdir();(root/"orphan"/"differential_expression.tsv").write_text("x\n");forge_seal(root);codes={q.code for q in self.release_report(root).diagnostics};self.assertIn("DUPLICATE_DATASET_ID",codes);self.assertIn("ORPHAN_DATASET_ARTIFACT",codes)
    def test_ghost_and_unauthorized_dataset_fail(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x));e=rows(root/"dataset_eligibility.tsv");ghost=dict(e[0]);ghost["dataset_id"]="GSE_FAKE";e.append(ghost);tsv(root/"dataset_eligibility.tsv",e);forge_seal(root);codes={q.code for q in self.release_report(root).diagnostics};self.assertIn("UNAUTHORIZED_ELIGIBILITY_DATASET",codes);self.assertIn("ELIGIBILITY_WITHOUT_DATASET_DIRECTORY",codes)
    def test_gse75665_de_is_valid_but_saturated_design_fails(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x),"GSE75665",2);self.assertEqual(self.dataset_report(d,"GSE75665").error_count,0)
            design=rows(d/"design.tsv")
            for row in design:row["formula"]="~ 0 + subject:altitude"
            tsv(d/"design.tsv",design)
            self.assertIn("GSE75665_SATURATED_DESIGN",{q.code for q in self.dataset_report(d,"GSE75665").diagnostics})
    def test_incomplete_selected_candidate_cannot_unlock_phase_b(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x),True);(root/"GSE103940"/"differential_expression.tsv").unlink();forge_seal(root);r=self.release_report(root);self.assertFalse(r.metrics["phase_b_complete"]);self.assertIn("MISSING_DATASET_ARTIFACT",{q.code for q in r.diagnostics})
    def test_meta_contributor_hash_and_effect_set_fail(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x),True);mm=json.loads((root/"meta_method.json").read_text());mm["contributors"][1]["source_de_sha256"]=H;(root/"meta_method.json").write_text(json.dumps(mm));ef=rows(root/"meta_study_effects.tsv");ef[0]["effect"]="999";ef=[r for r in ef if r["dataset_id"]!="GSE103940"];tsv(root/"meta_study_effects.tsv",ef);forge_seal(root);codes={q.code for q in self.release_report(root).diagnostics};self.assertIn("META_CONTRIBUTOR_HASH_MISMATCH",codes);self.assertIn("META_STUDY_EFFECT_COVERAGE_MISMATCH",codes);self.assertIn("META_STUDY_EFFECT_VALUE_MISMATCH",codes)
    def test_meta_multicontrast_source_family_selection(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x),True);d196=root/"GSE196728";add_second_total_bulk_family(d196);dirs={"GSE196728":d196,"GSE103940":root/"GSE103940"};summaries={}
            for ds,d in dirs.items():
                rr=V.Report("dataset");summaries[ds]=V.validate_dataset(d,ds,rr);self.assertEqual(rr.error_count,0,rr.payload())
            mm=json.loads((root/"meta_method.json").read_text());primary=mm["contributors"][0];primary["source_de_sha256"]=digest(d196/"differential_expression.tsv");primary["dataset_payload_sha256"]=summaries["GSE196728"]["dataset_payload_sha256"]
            effects=rows(root/"meta_study_effects.tsv")
            for r in effects:
                if r["dataset_id"]=="GSE196728":r["source_de_sha256"]=primary["source_de_sha256"]
            tsv(root/"meta_study_effects.tsv",effects);native=rows(root/"meta_native_results.tsv");published=rows(root/"altitude_whole_blood_meta.tsv")
            for target in (native,published):
                for r in target:r["source_effects_sha256"]=V.source_effects_sha256([e for e in effects if e["test_id"]==r["test_id"]])
            tsv(root/"meta_native_results.tsv",native);tsv(root/"altitude_whole_blood_meta.tsv",published);mm["native_results"]["sha256"]=digest(root/"meta_native_results.tsv");(root/"meta_method.json").write_text(json.dumps(mm))
            valid=V.Report("meta");V.validate_meta(root,dirs,summaries,["GSE196728","GSE103940"],valid);self.assertEqual(valid.error_count,0,valid.payload())
            mm=json.loads((root/"meta_method.json").read_text());mm["contributors"][0].update({"source_family_id":"fam_total_bulk_secondary","model_id":"m_total_bulk_secondary","contrast_id":"target_secondary","scientific_comparison_id":"altitude_secondary"});(root/"meta_method.json").write_text(json.dumps(mm));bad=V.Report("meta");V.validate_meta(root,dirs,summaries,["GSE196728","GSE103940"],bad);codes={q.code for q in bad.diagnostics};self.assertIn("META_CROSS_STUDY_COMPARISON_MISMATCH",codes);self.assertIn("META_STUDY_EFFECT_FOREIGN_KEY_MISMATCH",codes)
    def test_duplicate_design_keys_fail_in_all_orders_before_preseal(self):
        with tempfile.TemporaryDirectory() as x:
            base=Path(x)
            for comparison,prepend in (("altitude_primary",False),("altitude_primary",True),("contradictory_comparison",False),("contradictory_comparison",True)):
                with self.subTest(comparison=comparison,prepend=prepend):
                    root=build_release(base/f"{comparison}-{prepend}",False,False);d=root/"GSE196728";design=rows(d/"design.tsv");original=dict(design[0]);duplicate=dict(original);duplicate["scientific_comparison_id"]=comparison
                    design=([duplicate]+design) if prepend else (design+[duplicate]);tsv(d/"design.tsv",design);manifest=write_dataset_manifest(d,"GSE196728","WHOLE_BLOOD_DE")
                    eligibility=rows(root/"dataset_eligibility.tsv")
                    for row in eligibility:
                        if row["dataset_id"]=="GSE196728":
                            row["dataset_manifest_sha256"]=digest(d/"dataset_manifest.json");row["dataset_payload_sha256"]=manifest["dataset_payload_sha256"]
                    tsv(root/"dataset_eligibility.tsv",eligibility);build_index(root)
                    dataset=self.dataset_report(d,"GSE196728");self.assertIn("DUPLICATE_DESIGN_KEY",{q.code for q in dataset.diagnostics})
                    report=base/f"{comparison}-{prepend}.json";self.assertEqual(V.main(["release-preseal","--release-dir",str(root),"--report",str(report)]),1)
                    self.assertIn("DUPLICATE_DESIGN_KEY",{q["code"] for q in json.loads(report.read_text())["diagnostics"]});self.assertFalse((root/"validation_report.json").exists());self.assertFalse((root/"_SUCCESS.json").exists())
    def test_every_method_family_has_one_design_foreign_key(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x),"GSE196728",2);method=json.loads((d/"method.json").read_text());ghost=dict(method["test_families"][0]);ghost.update({"bh_family_id":"ghost_family","model_id":"ghost_model"});method["test_families"].append(ghost);(d/"method.json").write_text(json.dumps(method));write_dataset_manifest(d,"GSE196728","WHOLE_BLOOD_DE")
            self.assertIn("METHOD_FAMILY_DESIGN_FOREIGN_KEY_MISMATCH",{q.code for q in self.dataset_report(d,"GSE196728").diagnostics})
    def test_two_envelope_end_to_end_and_cross_payload_rejection(self):
        with tempfile.TemporaryDirectory() as x:
            base=Path(x);root=build_release(base/"release",False,False);pre=base/"preseal.json";final=base/"final.json";self.assertFalse((root/"validation_report.json").exists());self.assertFalse((root/"_SUCCESS.json").exists())
            self.assertEqual(V.main(["release-preseal","--release-dir",str(root),"--report",str(pre)]),0);submitted_bytes=pre.read_bytes();self.assertEqual(P.main(["--release-dir",str(root),"--preseal-report",str(pre)]),0);self.assertEqual(pre.read_bytes(),submitted_bytes);self.assertEqual(json.loads((root/"validation_report.json").read_text())["status"],"PASS");self.assertEqual(json.loads((root/"_SUCCESS.json").read_text())["validation_method"],"CONTROLLED_PRESEAL_RERUN");self.assertEqual(V.main(["release-final","--release-dir",str(root),"--report",str(final)]),0)
            bad=build_release(base/"bad",False,False);comp=rows(bad/"cell_composition_sensitivity.tsv");comp[0]["total_bulk_log2fc"]="999";tsv(bad/"cell_composition_sensitivity.tsv",comp);build_index(bad);bad_report=base/"bad-preseal.json";self.assertEqual(V.main(["release-preseal","--release-dir",str(bad),"--report",str(bad_report)]),1);self.assertEqual(json.loads(bad_report.read_text())["status"],"FAIL");self.assertEqual(P.main(["--release-dir",str(bad),"--preseal-report",str(bad_report)]),1);self.assertFalse((bad/"_SUCCESS.json").exists())
            other=build_release(base/"other",True,False);self.assertEqual(P.main(["--release-dir",str(other),"--preseal-report",str(pre)]),1);shutil.copyfile(pre,other/"validation_report.json");other_payload=hashlib.sha256(V.canonical_index(rows(other/"artifact_index.tsv")).encode()).hexdigest();marker(other,other_payload);cross=V.Report("release");V.validate_release(other,cross,"final");self.assertIn("INVALID_SEALED_VALIDATION_REPORT",{q.code for q in cross.diagnostics})
    def test_sealer_rejects_every_report_contract_mutation_and_fail_rewrite(self):
        with tempfile.TemporaryDirectory() as x:
            base=Path(x)
            mutations={
                "error_count":lambda p:p.__setitem__("error_count",7),
                "diagnostics":lambda p:p["diagnostics"].append({"severity":"ERROR","code":"FORGED","message":"forged","path":"","row":None}),
                "metrics":lambda p:p["metrics"].__setitem__("phase_b_complete",not p["metrics"]["phase_b_complete"]),
                "phase":lambda p:p.__setitem__("phase_b_complete",not p["phase_b_complete"]),
                "publishable":lambda p:p.__setitem__("publishable",True),
                "time":lambda p:(p.__setitem__("generated_at_utc","2000-01-01T00:00:00+00:00"),p["generation_context"].__setitem__("generated_at_utc","2000-01-01T00:00:00+00:00")),
                "command":lambda p:p["generation_context"].__setitem__("command","release-preseal --forged"),
                "extra":lambda p:p.__setitem__("forged_extra",True),
            }
            for name,mutate in mutations.items():
                with self.subTest(name=name):
                    root=build_release(base/name,False,False);pre=base/f"{name}.json";self.assertEqual(V.main(["release-preseal","--release-dir",str(root),"--report",str(pre)]),0);payload=json.loads(pre.read_text());mutate(payload);pre.write_bytes(V.report_bytes(payload))
                    self.assertEqual(P.main(["--release-dir",str(root),"--preseal-report",str(pre)]),1);self.assertFalse((root/"validation_report.json").exists());self.assertFalse((root/"_SUCCESS.json").exists())
            bad=build_release(base/"fail-rewrite",False,False);composition=rows(bad/"cell_composition_sensitivity.tsv");composition[0]["total_bulk_log2fc"]="999";tsv(bad/"cell_composition_sensitivity.tsv",composition);build_index(bad);failed=base/"failed.json";self.assertEqual(V.main(["release-preseal","--release-dir",str(bad),"--report",str(failed)]),1);payload=json.loads(failed.read_text());payload.update({"status":"PASS","schema_pass":True,"scientific_gate_pass":True});failed.write_bytes(V.report_bytes(payload))
            self.assertEqual(P.main(["--release-dir",str(bad),"--preseal-report",str(failed)]),1);self.assertFalse((bad/"validation_report.json").exists());self.assertFalse((bad/"_SUCCESS.json").exists())
    def test_sealer_rolls_back_partial_envelope_on_commit_failure(self):
        with tempfile.TemporaryDirectory() as x:
            base=Path(x);root=build_release(base/"release",False,False);pre=base/"preseal.json";self.assertEqual(V.main(["release-preseal","--release-dir",str(root),"--report",str(pre)]),0);real_replace=P.os.replace;calls=0
            def fail_second(src,dst):
                nonlocal calls
                calls+=1
                if calls==1:return real_replace(src,dst)
                raise OSError("simulated marker commit failure")
            with mock.patch.object(P.os,"replace",side_effect=fail_second):
                self.assertEqual(P.main(["--release-dir",str(root),"--preseal-report",str(pre)]),1)
            self.assertFalse((root/"validation_report.json").exists());self.assertFalse((root/"_SUCCESS.json").exists())
    def test_final_reconstructs_full_report_and_exact_marker_contracts(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x));report_path=root/"validation_report.json";marker_path=root/"_SUCCESS.json";original_report=report_path.read_bytes();original_marker=marker_path.read_bytes()
            mutations={
                "diagnostics":lambda p:p["diagnostics"].append({"severity":"ERROR","code":"FORGED","message":"forged","path":"","row":None}),
                "metrics":lambda p:p["metrics"].__setitem__("phase_b_complete",not p["metrics"]["phase_b_complete"]),
                "phase":lambda p:p.__setitem__("phase_b_complete",not p["phase_b_complete"]),
                "publishable":lambda p:p.__setitem__("publishable",True),
                "time":lambda p:p.__setitem__("generated_at_utc","2000-01-01T00:00:00+00:00"),
                "command":lambda p:p["generation_context"].__setitem__("command","release-preseal --forged"),
                "extra":lambda p:p.__setitem__("forged_extra",True),
            }
            for name,mutate in mutations.items():
                with self.subTest(name=name):
                    report_path.write_bytes(original_report);marker_path.write_bytes(original_marker);payload=json.loads(original_report);mutate(payload);report_path.write_bytes(V.report_bytes(payload));success=json.loads(original_marker);success["validation_report_sha256"]=digest(report_path);marker_path.write_bytes(V.report_bytes(success))
                    checked=self.release_report(root);self.assertIn("INVALID_SEALED_VALIDATION_REPORT",{q.code for q in checked.diagnostics})
            report_path.write_bytes(original_report);success=json.loads(original_marker);success["forged_extra"]=True;marker_path.write_bytes(V.report_bytes(success));self.assertIn("SUCCESS_MARKER_MISMATCH",{q.code for q in self.release_report(root).diagnostics})
    def test_renv_standard_dependency_fields_require_closure(self):
        with tempfile.TemporaryDirectory() as x:
            d=build_dataset(Path(x));lock=json.loads((d/"renv.lock").read_text());del lock["Packages"]["BiocGenerics"];(d/"renv.lock").write_text(json.dumps(lock));method=json.loads((d/"method.json").read_text());method["environment_lock"]["sha256"]=digest(d/"renv.lock");(d/"method.json").write_text(json.dumps(method));self.assertIn("LOCK_DEPENDENCY_MISSING",{q.code for q in self.dataset_report(d,"GSE46480").diagnostics})
    def test_composition_values_and_classification_are_source_derived(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x),True);comp=rows(root/"cell_composition_sensitivity.tsv");comp[0].update({"total_bulk_log2fc":"999","conditional_log2fc":"-999","delta_log2fc":"0","sensitivity_class":"direction_changed"});tsv(root/"cell_composition_sensitivity.tsv",comp);forge_seal(root);codes={q.code for q in self.release_report(root).diagnostics};self.assertIn("COMPOSITION_VALUE_MISMATCH",codes);self.assertIn("COMPOSITION_CLASS_MISMATCH",codes)
    def test_pbmc_target_provenance_and_statistics_are_bound(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x));pb=rows(root/"altitude_pbmc_validation.tsv");pb[0].update({"target_source_sha256":H,"statistic":"999"});tsv(root/"altitude_pbmc_validation.tsv",pb);forge_seal(root);codes={q.code for q in self.release_report(root).diagnostics};self.assertIn("PBMC_TARGET_FOREIGN_KEY_MISMATCH",codes);self.assertIn("PBMC_VALUE_MISMATCH",codes)
    def test_meta_extreme_numbers_and_sensitivities_fail(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x),True);meta=rows(root/"altitude_whole_blood_meta.tsv");meta[0].update({"meta_log2fc":"999","meta_se":"0.000001","p_value":"0.9","q_heterogeneity":"999","tau2":"999","modified_hksj_low":"999"});tsv(root/"altitude_whole_blood_meta.tsv",meta);forge_seal(root);codes={q.code for q in self.release_report(root).diagnostics};self.assertIn("META_RESULT_VALUE_MISMATCH",codes);self.assertIn("META_FDR_MISMATCH",codes)
    def test_forged_sealed_report_fails_even_when_resealed(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x));sealed=json.loads((root/"validation_report.json").read_text());sealed["scientific_gate_pass"]=False;(root/"validation_report.json").write_text(json.dumps(sealed));success=json.loads((root/"_SUCCESS.json").read_text());success["validation_report_sha256"]=digest(root/"validation_report.json");(root/"_SUCCESS.json").write_text(json.dumps(success));self.assertIn("INVALID_SEALED_VALIDATION_REPORT",{q.code for q in self.release_report(root).diagnostics})
    def test_report_inside_payload_fails_without_write(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x));out=root/"inside.json";self.assertEqual(V.main(["release","--release-dir",str(root),"--report",str(out)]),1);self.assertFalse(out.exists())
    def test_cli_external_report_is_read_only(self):
        with tempfile.TemporaryDirectory() as x:
            root=build_release(Path(x)/"release");before={p.relative_to(root).as_posix():digest(p) for p in root.rglob("*") if p.is_file()};out=Path(x)/"out.json";self.assertEqual(V.main(["release","--release-dir",str(root),"--report",str(out)]),0);self.assertEqual(before,{p.relative_to(root).as_posix():digest(p) for p in root.rglob("*") if p.is_file()})
    def test_medication_closed_mapping(self):
        self.assertEqual([V.medication01(x) for x in (False,0,"no",True,1,"YES")],[0,0,0,1,1,1]);self.assertRaises(ValueError,V.medication01,"unknown")
    def test_r_and_python_medication_contract_are_synchronized(self):
        reference=SCRIPT.parents[2]/".claude/skills/module2-geo-expression/references/design-and-de.md"
        if not reference.is_file(): self.skipTest("internal design reference not present in this tree")
        text=reference.read_text()
        for token in ('"false" = 0','"0" = 0','"no" = 0','"true" = 1','"1" = 1','"yes" = 1'):self.assertIn(token,text)
        self.assertNotIn('allowed <- c("no", "yes")',text)

if __name__=="__main__":unittest.main()
