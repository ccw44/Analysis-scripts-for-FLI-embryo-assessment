#!/usr/bin/env python
"""Supplementary checks: published loci re-checked on mouse mm39 coordinates, coverage
sensitivity of the analysable set, and agreement between the DSS and region-level calls.
"""
import os
import numpy as np
import pandas as pd

# Data locations - edit this block only; run from the repository root.
#   IN_CPG_DIR  per-sample CpG tables: chrom, start0, end, beta, M, U (TSV, no header)
#   REF_DIR     mm39: mm39.chrom.sizes, refGene.txt.gz, cpgIslandExt.txt.gz
#   OUT_DIR     output root: tables/ and figures/
IN_CPG_DIR = "data/cpg_tables"
REF_DIR = "data/mm39"
OUT_DIR = "results"
TAB = os.path.join(OUT_DIR, "tables")
FIG = os.path.join(OUT_DIR, "figures")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
REF = REF_DIR

ORDER = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"]
PROM_UP, PROM_DOWN = 2000, 500

# ---------- (a) published loci on mouse mm39 ---------------------------------
cols = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
        "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat",
        "cdsEndStat", "exonFrames"]
rg = pd.read_csv(os.path.join(REF, "refGene.txt.gz"), sep="\t", header=None, names=cols,
                 low_memory=False)
TARGETS = ["Dnmt1", "Dnmt3a", "Dnmt3b", "Tet1", "Tet2", "Tet3"]
rows = []
for tissue, ctrl, far in (
    ("liver", ["Control-liver-1", "Control-liver-2", "Control-liver-3"],
              ["FAR-liver-1", "FAR-liver-2", "FAR-liver-3"]),
    ("blood", ["Control-Blood-1", "Control-Blood-2", "Control-Blood-3"],
              ["FAR-Blood-1", "FAR-Blood-2", "FAR-Blood-3"]),
):
    comp = pd.read_csv(os.path.join(TAB, f"comparable_{tissue}.tsv.gz"), sep="\t")
    ann = pd.read_csv(os.path.join(TAB, f"DMR_annotation_{tissue}.tsv.gz"), sep="\t",
                      usecols=["region_id", "chrom", "start", "end", "q_primary", "diff"])
    for sym in TARGETS:
        sub = rg[rg["name2"] == sym]
        if len(sub) == 0:
            continue
        for chrom in sorted(sub["chrom"].unique(), key=lambda c: ORDER.index(c)):
            s2 = sub[sub["chrom"] == chrom]
            wins = []
            for _, r in s2.iterrows():
                if r["strand"] == "+":
                    wins.append((r["txStart"] - PROM_UP, r["txStart"] + PROM_DOWN))
                else:
                    wins.append((r["txEnd"] - PROM_DOWN, r["txEnd"] + PROM_UP))
            for ws, we in sorted(set(wins)):
                m = (comp["chrom"] == chrom) & (comp["start"] >= ws) & (comp["end"] <= we)
                if m.sum() == 0:
                    continue
                sel = comp[m]
                vals = {}
                for s in ctrl + far:
                    b = (sel[s + "_M"] / sel[s + "_cov"]).replace([np.inf, -np.inf], np.nan)
                    vals[s] = float(np.nanmean(b)) if b.notna().any() else np.nan
                cm = np.nanmean([vals[s] for s in ctrl])
                fm = np.nanmean([vals[s] for s in far])
                ov = ann[(ann["chrom"] == chrom) & (ann["start"] < we) & (ann["end"] > ws)]
                rows.append(dict(
                    tissue=tissue, gene=sym, chrom=chrom,
                    window_start=int(ws), window_end=int(we),
                    n_comparable_cpg=int(m.sum()),
                    n_transcripts=int(len(s2)),
                    control_mean=cm, far_mean=fm, diff=fm - cm,
                    n_candidate_regions_overlapping=int(len(ov)),
                    min_region_q=float(ov["q_primary"].min()) if len(ov) else np.nan))
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(TAB, "supplementary_loci_check.csv"), index=False, float_format="%.6g")
print("=== published loci (mouse mm39 promoter windows, CpGs from the comparable set) ===")
print(tab.to_string(index=False, max_colwidth=24))

# ---------- (b) coverage sensitivity ----------------------------------------
sens = []
for tissue, ctrl, far in (
    ("liver", ["Control-liver-1", "Control-liver-2", "Control-liver-3"],
              ["FAR-liver-1", "FAR-liver-2", "FAR-liver-3"]),
    ("blood", ["Control-Blood-1", "Control-Blood-2", "Control-Blood-3"],
              ["FAR-Blood-1", "FAR-Blood-2", "FAR-Blood-3"]),
):
    comp = pd.read_csv(os.path.join(TAB, f"comparable_{tissue}.tsv.gz"), sep="\t",
                       usecols=["chrom", "start"] + [s + "_cov" for s in ctrl + far])
    ctrl_cnt = (comp[[s + "_cov" for s in ctrl]] >= 10).sum(axis=1)
    far_cnt = (comp[[s + "_cov" for s in far]] >= 10).sum(axis=1)
    m10 = (ctrl_cnt >= 2) & (far_cnt >= 2)
    pos = comp.loc[m10, "start"].to_numpy()
    chr_ = comp.loc[m10, "chrom"].to_numpy()
    gaps = []
    for c in np.unique(chr_):
        p = np.sort(pos[chr_ == c])
        if len(p) > 1:
            gaps.append(np.diff(p))
    g = np.concatenate(gaps) if gaps else np.array([np.nan])
    sens.append(dict(tissue=tissue,
                     n_comparable_ge5=len(comp), n_comparable_ge10=int(m10.sum()),
                     frac_of_ge5=round(float(m10.mean()), 5),
                     median_spacing_ge5_bp=float(np.median(np.diff(np.sort(pos)))) if m10.sum() > 1 else np.nan,
                     median_spacing_ge10_bp=float(np.median(g)),
                     frac_gap_lt_100bp_ge10=round(float(np.mean(g < 100)), 4)))
sens = pd.DataFrame(sens)
sens.to_csv(os.path.join(TAB, "coverage_sensitivity.csv"), index=False, float_format="%.4g")
print("\n=== coverage sensitivity ===")
print(sens.to_string(index=False))

# ---------- (c) DSS vs region-level agreement -------------------------------
agree = []
for tissue in ("liver", "blood"):
    dss = pd.read_csv(os.path.join(TAB, f"DSS_regions_{tissue}.tsv"), sep="\t")
    ann = pd.read_csv(os.path.join(TAB, f"DMR_annotation_{tissue}.tsv.gz"), sep="\t",
                      usecols=["region_id", "chrom", "start", "end", "q_primary", "diff"])
    for _, r in dss.iterrows():
        m = (ann["chrom"] == r["chr"]) & (ann["start"] < r["end"]) & (ann["end"] > r["start"])
        ov = ann[m]
        agree.append(dict(
            tissue=tissue, dss_chr=r["chr"], dss_start=int(r["start"]), dss_end=int(r["end"]),
            dss_nCG=int(r["nCG"]), dss_diff=float(r["diff.Methy"]),
            dss_areaStat=float(r["areaStat"]),
            n_overlapping_candidate_regions=int(len(ov)),
            min_region_q=float(ov["q_primary"].min()) if len(ov) else np.nan,
            max_abs_region_diff=float(ov["diff"].abs().max()) if len(ov) else np.nan,
            has_region_q05=bool((ov["q_primary"] < 0.05).any()) if len(ov) else False))
ag = pd.DataFrame(agree)
ag.to_csv(os.path.join(TAB, "dss_vs_region_agreement.csv"), index=False, float_format="%.6g")
print("\n=== DSS regions vs region-level inference ===")
print(ag.to_string(index=False, max_colwidth=18))
