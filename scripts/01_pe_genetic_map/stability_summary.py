#!/usr/bin/env python3
"""Frozen Module 1 stability-summary statistical domains and computation."""
from collections import defaultdict
import math

NA="NA"
DOMAINS={
    ("whole_blood_gene","signed_z"):("primary_stat","comparison_stat","signed"),
    ("whole_blood_gene","abs_z"):("primary_stat","comparison_stat","absolute"),
    ("whole_blood_candidate_projection","candidate_projection"):None,
    ("tissue_candidate_projection","candidate_projection"):None,
    ("locus","locus_direction"):None,
}

def _ranks(values):
    order=sorted(range(len(values)),key=lambda i:values[i]); result=[0.]*len(values); i=0
    while i<len(values):
        j=i+1
        while j<len(values) and values[order[j]]==values[order[i]]: j+=1
        rank=(i+j-1)/2+1
        for k in range(i,j): result[order[k]]=rank
        i=j
    return result

def _spearman(x,y):
    rx,ry=_ranks(x),_ranks(y); mx=sum(rx)/len(rx); my=sum(ry)/len(ry)
    vx=sum((a-mx)**2 for a in rx); vy=sum((b-my)**2 for b in ry)
    if vx==0 or vy==0: return None,"undefined_constant_rank"
    return sum((a-mx)*(b-my) for a,b in zip(rx,ry))/math.sqrt(vx*vy),"success"

def summarize(rows):
    grouped=defaultdict(list)
    for row in rows: grouped[(row["comparison_dataset"],row["comparison_level"],row["analysis_set"],row["rank_metric"])].append(row)
    output=[]
    for (dataset,level,analysis_set,metric), group in sorted(grouped.items()):
        domain=DOMAINS.get((level,metric))
        if (level,metric) not in DOMAINS: raise ValueError(f"unsupported stability summary domain: {level}/{metric}")
        numeric=[]
        if domain:
            left,right,transform=domain
            for row in group:
                if row[left]!=NA and row[right]!=NA:
                    pair=(float(row[left]),float(row[right])); numeric.append(tuple(abs(x) for x in pair) if transform=="absolute" else pair)
            if len(numeric)<10: rho,status=None,"insufficient_n"
            else: rho,status=_spearman(*map(list,zip(*numeric)))
        else: rho,status=None,"not_applicable_no_numeric_rank"
        output.append({"comparison_dataset":dataset,"comparison_level":level,"analysis_set":analysis_set,"rank_metric":metric,"n_common":len(numeric),"spearman_rho":rho,"status":status})
    return output

def discrepancies(observed, expected, tolerance=1e-12):
    key=lambda r:(r["comparison_dataset"],r["comparison_level"],r["analysis_set"],r["rank_metric"])
    obs={key(r):r for r in observed}; exp={key(r):r for r in expected}; errors=[]
    if len(obs)!=len(observed) or set(obs)!=set(exp): errors.append("summary key closure")
    for k in set(obs)&set(exp):
        o,e=obs[k],exp[k]
        if str(o["n_common"])!=str(e["n_common"]): errors.append(f"{k}: n_common")
        if o["status"]!=e["status"]: errors.append(f"{k}: status")
        expected_rho=e["spearman_rho"]
        if expected_rho is None:
            if o["spearman_rho"]!=NA: errors.append(f"{k}: rho")
        else:
            try:
                if not math.isfinite(float(o["spearman_rho"])) or abs(float(o["spearman_rho"])-expected_rho)>tolerance: errors.append(f"{k}: rho")
            except ValueError: errors.append(f"{k}: rho")
    return errors
