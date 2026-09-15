#!/usr/bin/env python3

"""GTEx v11 Whole Blood sQTL/apaQTL audit for WARS1 (post-selection, 2026-08-02).

Frozen outputs live in
results/04_core_candidates/analyses/m4_wars1_post_selection_audit_v1_20260802/isoform_audit/.
Reruns write to work/ and must be COMPARED against the frozen outputs, never
written back into the release.
"""

from datetime import date
from pathlib import Path
import os

import pandas as pd
import pyarrow.parquet as pq


OUT = Path("work/04_core_candidates/reruns/wars1_isoform_qtl") / date.today().strftime("%Y%m%d")
OUT.mkdir(parents=True, exist_ok=True)
# Local data root containing GTEx/ subdirectories downloaded from the public
# sources; override with the PE_DATA_ROOT environment variable.
DATA_ROOT = Path(os.environ.get("PE_DATA_ROOT", "data"))
SQTL = DATA_ROOT / "GTEx/GTEx_Analysis_v11_sQTL/Whole_Blood.v11.sQTLs.signif_pairs.parquet"
APA = DATA_ROOT / "GTEx/GTEx_Analysis_v11_apaQTL/Whole_Blood.v11.apaQTLs.signif_pairs.parquet"
GROUP = "ENSG00000140105.18"
LEAD_EQTL = "chr14_100374476_T_C_b38"


sqtl = pq.read_table(SQTL, filters=[("group_id", "==", GROUP)]).to_pandas()
apa = pq.read_table(APA, filters=[("group_id", "==", GROUP)]).to_pandas()

phenotype_summary = (
    sqtl.groupby("phenotype_id", as_index=False)
    .agg(n_significant_pairs=("variant_id", "size"), min_p=("pval_nominal", "min"))
    .sort_values("min_p")
)
lead = sqtl.loc[
    sqtl["variant_id"].eq(LEAD_EQTL),
    ["phenotype_id", "variant_id", "pval_nominal", "slope", "slope_se"],
].sort_values("pval_nominal")

phenotype_summary.to_csv(OUT / "wars1_whole_blood_v11_sqtl_phenotypes.tsv", sep="\t", index=False)
lead.to_csv(OUT / "wars1_lead_eqtl_as_sqtl.tsv", sep="\t", index=False)

overview = pd.DataFrame(
    {
        "metric": [
            "significant_sQTL_pairs",
            "sQTL_intron_phenotypes",
            "lead_eQTL_significant_sQTL_phenotypes",
            "significant_apaQTL_pairs",
        ],
        "value": [len(sqtl), sqtl["phenotype_id"].nunique(), len(lead), len(apa)],
    }
)
overview.to_csv(OUT / "wars1_isoform_qtl_overview.tsv", sep="\t", index=False)
print(overview.to_string(index=False))
print(lead.to_string(index=False))
