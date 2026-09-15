#!/usr/bin/env python3
"""Seal one already validated Module 2 payload; this does not publish or rename it.

⚠ LEGACY / DO NOT USE FOR CURRENT DELIVERABLES (2026-09-10, audit B-10): this
sealer wraps validate_module2.py, which implements the retired 5.0 payload
schema and a dataset whitelist of four now-excluded datasets. The current
GSE103927/GSE333506 v2 deliverables use the per-analysis freeze discipline
(freeze_registry.json + source_manifest.tsv sha256 closure) and
validate_module2_v2_bh.py instead.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
from datetime import datetime
from pathlib import Path

VALIDATOR_PATH=Path(__file__).with_name("validate_module2.py")
SPEC=importlib.util.spec_from_file_location("module2_validator_for_seal",VALIDATOR_PATH)
assert SPEC and SPEC.loader
V=importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name]=V
SPEC.loader.exec_module(V)

def seal(release_dir:Path,preseal_report:Path):
    root=release_dir.resolve();source=preseal_report.resolve();sealed=root/"validation_report.json";marker=root/"_SUCCESS.json"
    if sealed.exists() or marker.exists():raise ValueError("refusing to overwrite an existing envelope")
    try:source.relative_to(root)
    except ValueError:pass
    else:raise ValueError("pre-seal report must be outside the payload")
    source_bytes=source.read_bytes();pre=json.loads(source_bytes)
    if not isinstance(pre,dict) or not V.valid_utc(pre.get("generated_at_utc")):raise ValueError("pre-seal report has no valid generation time")
    generated_at=str(pre["generated_at_utc"])
    generated_epoch=datetime.fromisoformat(generated_at.replace("Z","+00:00")).timestamp()
    if abs(source.stat().st_mtime-generated_epoch)>5:raise ValueError("pre-seal report timestamp does not match the report file")
    submitted_replay=V.controlled_preseal_report(root,generated_at)
    if submitted_replay.error_count or source_bytes!=V.report_bytes(submitted_replay.payload()):
        raise ValueError("report is not the complete current-validator pre-seal output")

    authoritative=V.controlled_preseal_report(root)
    if authoritative.error_count:raise ValueError("controlled pre-seal validator rerun failed")
    authoritative_bytes=V.report_bytes(authoritative.payload())
    index_check=V.Report("publisher-index-check");index_hash,payload_hash=V.validate_artifact_index(root,index_check)
    if index_check.error_count or authoritative.metrics.get("payload_root_sha256")!=payload_hash:
        raise ValueError("payload changed during controlled validation")
    success={"payload_root_sha256":payload_hash,"artifact_index_sha256":index_hash,"validator_sha256":V.sha(VALIDATOR_PATH),"sealer_sha256":V.sha(Path(__file__)),"schema_version":V.SCHEMA_VERSION,"schema_sha256":V.schema_sha256(),"validation_report_sha256":hashlib.sha256(authoritative_bytes).hexdigest(),"validation_method":"CONTROLLED_PRESEAL_RERUN"}
    success_bytes=(json.dumps(success,indent=2,sort_keys=True)+"\n").encode("utf-8")
    committed=False
    try:
        with tempfile.TemporaryDirectory(prefix=".module2-seal-",dir=root.parent) as tmp:
            tmp_root=Path(tmp);report_tmp=tmp_root/"validation_report.json";marker_tmp=tmp_root/"_SUCCESS.json"
            report_tmp.write_bytes(authoritative_bytes);marker_tmp.write_bytes(success_bytes)
            os.replace(report_tmp,sealed)
            os.replace(marker_tmp,marker)
            committed=True
    finally:
        if not committed:
            marker.unlink(missing_ok=True)
            sealed.unlink(missing_ok=True)

def main(argv=None):
    parser=argparse.ArgumentParser();parser.add_argument("--release-dir",type=Path,required=True);parser.add_argument("--preseal-report",type=Path,required=True);args=parser.parse_args(argv)
    try:seal(args.release_dir,args.preseal_report)
    except (OSError,ValueError,json.JSONDecodeError) as exc:
        print(f"SEAL FAIL: {exc}",file=sys.stderr);return 1
    print("SEAL PASS");return 0

if __name__=="__main__":sys.exit(main())
