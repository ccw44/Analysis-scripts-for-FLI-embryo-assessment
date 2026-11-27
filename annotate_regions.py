#!/usr/bin/env python
"""Region annotation, result layers and per-tissue summary counts.

Feature classes are mutually exclusive (promoter > exon > intron > intergenic); CpG-island
classes use a separate exclusive scheme (island > shore > shelf > open sea). Result layers
are candidate regions, q < 0.05 and the large-effect subset; an empty layer is still
written, with a header.
"""
import os, gzip
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
CIDX = {c: i for i, c in enumerate(ORDER)}
PROM_UP, PROM_DOWN = 2000, 500
SHORE, SHELF = 2000, 4000

# ---------------- reference tables ----------------
def load_refgene():
    cols = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
            "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat",
            "cdsEndStat", "exonFrames"]
    df = pd.read_csv(os.path.join(REF, "refGene.txt.gz"), sep="\t", header=None, names=cols,
                     dtype={"chrom": "string"}, low_memory=False)
    df = df[df["chrom"].isin(ORDER)].copy()
    return df

def load_cgi():
    df = pd.read_csv(os.path.join(REF, "cpgIslandExt.txt.gz"), sep="\t", header=None,
                     names=["bin", "chrom", "chromStart", "chromEnd", "name", "length",
                            "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"],
                     dtype={"chrom": "string"})
    return df[df["chrom"].isin(ORDER)].copy()

def merge_intervals(starts, ends):
    """merge into sorted, non-overlapping intervals"""
    o = np.argsort(starts, kind="stable")
    s = np.asarray(starts)[o]; e = np.asarray(ends)[o]
    ms, me = [], []
    cs, ce = s[0], e[0]
    for a, b in zip(s[1:], e[1:]):
        if a <= ce:
            ce = max(ce, b)
        else:
            ms.append(cs); me.append(ce); cs, ce = a, b
    ms.append(cs); me.append(ce)
    return np.array(ms), np.array(me)

def overlap_mask(qs, qe, ms, me):
    """qs/qe query intervals; ms/me merged non-overlapping reference intervals"""
    if len(ms) == 0:
        return np.zeros(len(qs), dtype=bool)
    i = np.searchsorted(ms, qe, side="left") - 1        # last interval starting before qe
    ok = np.zeros(len(qs), dtype=bool)
    valid = i >= 0
    ok[valid] = me[i[valid]] > qs[valid]
    return ok

def distance_to_nearest(qs, qe, ts):
    """signed distance from query midpoint to nearest TSS (ts sorted)"""
    if len(ts) == 0:
        return np.full(len(qs), np.nan), np.zeros(len(qs), dtype=np.int64)
    mid = (qs + qe) // 2
    i = np.searchsorted(ts, mid)
    left = np.clip(i - 1, 0, len(ts) - 1)
    right = np.clip(i, 0, len(ts) - 1)
    dl = np.abs(mid - ts[left]); dr = np.abs(mid - ts[right])
    use_l = dl <= dr
    dist = np.where(use_l, mid - ts[left], mid - ts[right])
    j = np.where(use_l, left, right)
    return dist, j

rg = load_refgene()
cgi = load_cgi()
print(f"refGene transcripts: {len(rg):,}  genes: {rg['name2'].nunique():,}")
print(f"CpG islands: {len(cgi):,}")

# per-chromosome reference structures
EXON, PROM, BODY, TSS, TSSGENE, CGI, GENES = {}, {}, {}, {}, {}, {}, {}
for c in ORDER:
    sub = rg[rg["chrom"] == c]
    if len(sub) == 0:
        continue
    ex_s, ex_e = [], []
    for s, e in zip(sub["exonStarts"], sub["exonEnds"]):
        ss = np.fromstring(s.rstrip(","), sep=",", dtype=np.int64)
        ee = np.fromstring(e.rstrip(","), sep=",", dtype=np.int64)
        ex_s.append(ss); ex_e.append(ee)
    ex_s = np.concatenate(ex_s); ex_e = np.concatenate(ex_e)
    EXON[c] = merge_intervals(ex_s, ex_e)
    BODY[c] = merge_intervals(sub["txStart"].to_numpy(), sub["txEnd"].to_numpy())
    tss = np.where(sub["strand"].to_numpy() == "+", sub["txStart"].to_numpy(), sub["txEnd"].to_numpy())
    TSS[c] = np.sort(tss)
    tssg = np.where(sub["strand"].to_numpy() == "+", sub["txStart"].to_numpy(), sub["txEnd"].to_numpy())
    o = np.argsort(tssg, kind="stable")
    TSSGENE[c] = sub["name2"].to_numpy()[o]
    ps = np.where(sub["strand"].to_numpy() == "+", sub["txStart"].to_numpy() - PROM_UP,
                  sub["txEnd"].to_numpy() - PROM_DOWN)
    pe = np.where(sub["strand"].to_numpy() == "+", sub["txStart"].to_numpy() + PROM_DOWN,
                  sub["txEnd"].to_numpy() + PROM_UP)
    PROM[c] = merge_intervals(ps, pe)
    # gene bodies with names for multi-gene overlap
    gs = sub["txStart"].to_numpy(); ge = sub["txEnd"].to_numpy()
    go = np.argsort(gs, kind="stable")
    GENES[c] = (gs[go], ge[go], sub["name2"].to_numpy()[go],
                np.maximum.accumulate(ge[go]))
    cs = cgi[cgi["chrom"] == c]
    CGI[c] = merge_intervals(cs["chromStart"].to_numpy(), cs["chromEnd"].to_numpy()) if len(cs) else (np.array([]), np.array([]))

def genes_overlapping(c, qs, qe):
    """names of all genes whose body overlaps each query (O(log n + k))"""
    if c not in GENES:
        return [[] for _ in range(len(qs))]
    gs, ge, gn, gmax = GENES[c]
    out = []
    for s, e in zip(qs, qe):
        i = np.searchsorted(gs, e, side="left") - 1
        names = []
        while i >= 0 and gmax[i] > s:
            if ge[i] > s:
                names.append(gn[i])
            i -= 1
        out.append(names)
    return out

def classify(chrom_arr, qs, qe):
    n = len(qs)
    feat = np.empty(n, dtype=object)
    cgi_cls = np.empty(n, dtype=object)
    nearest_gene = np.empty(n, dtype=object)
    nearest_dist = np.full(n, np.nan)
    n_genes = np.zeros(n, dtype=np.int32)
    gene_list = np.empty(n, dtype=object)
    for c in pd.unique(chrom_arr):
        m = np.flatnonzero(chrom_arr == c)
        s, e = qs[m], qe[m]
        in_prom = overlap_mask(s, e, *PROM.get(c, (np.array([]), np.array([]))))
        in_exon = overlap_mask(s, e, *EXON.get(c, (np.array([]), np.array([]))))
        in_body = overlap_mask(s, e, *BODY.get(c, (np.array([]), np.array([]))))
        cls = np.where(in_prom, "promoter", np.where(in_exon, "exon",
              np.where(in_body, "intron", "intergenic")))
        feat[m] = cls
        # CpG island class
        isl = CGI.get(c, (np.array([]), np.array([])))
        in_isl = overlap_mask(s, e, *isl)
        if len(isl[0]):
            d_is = np.full(len(s), np.inf)
            mid = (s + e) // 2
            i = np.searchsorted(isl[0], mid)
            lo = np.clip(i - 1, 0, len(isl[0]) - 1); hi = np.clip(i, 0, len(isl[0]) - 1)
            dl = np.where(mid < isl[0][lo], isl[0][lo] - mid, np.where(mid > isl[1][lo], mid - isl[1][lo], 0))
            dh = np.where(mid < isl[0][hi], isl[0][hi] - mid, np.where(mid > isl[1][hi], mid - isl[1][hi], 0))
            d_is = np.minimum(dl, dh)
            cc = np.where(in_isl, "island",
                 np.where(d_is <= SHORE, "shore",
                 np.where(d_is <= SHELF, "shelf", "open_sea")))
        else:
            cc = np.array(["open_sea"] * len(s), dtype=object)
        cgi_cls[m] = cc
        if c in TSS:
            dist, j = distance_to_nearest(s, e, TSS[c])
            nearest_dist[m] = dist
            nearest_gene[m] = TSSGENE[c][j]
        else:
            nearest_gene[m] = None
        gl = genes_overlapping(c, s, e)
        n_genes[m] = [len(x) for x in gl]
        gene_list[m] = [";".join(sorted(set(x))) if x else "" for x in gl]
    return feat, cgi_cls, nearest_gene, nearest_dist, n_genes, gene_list

summary_rows, dist_rows = [], []
for tissue in ("liver", "blood"):
    print(f"\n=== {tissue} ===", flush=True)
    reg = pd.read_csv(os.path.join(TAB, f"regions_all_{tissue}.tsv.gz"), sep="\t")
    test = pd.read_csv(os.path.join(TAB, f"region_tests_{tissue}.tsv.gz"), sep="\t")
    df = reg.merge(test, on="region_id", how="left", validate="one_to_one")
    print(f"regions={len(df):,}  primary q<0.05={int((df['q_primary']<0.05).sum())}  "
          f"sensitivity q<0.05={int((df['q_sensitivity']<0.05).sum())}")

    chrom_arr = df["chrom"].to_numpy()
    qs = df["start"].to_numpy(); qe = df["end"].to_numpy()
    feat, cgi_cls, ng, nd, ngene, glist = classify(chrom_arr, qs, qe)
    df["feature_class"] = feat
    df["cpg_island_class"] = cgi_cls
    df["nearest_gene"] = ng
    df["distance_to_TSS"] = nd          # genomic midpoint minus TSS, + = higher coordinate
    df["n_overlapping_genes"] = ngene
    df["overlapping_genes"] = glist
    cov_cols = [c for c in df.columns if c.endswith("_mean_cov")]
    df["mean_cov_region"] = df[cov_cols].mean(axis=1) if cov_cols else np.nan

    # result layers
    df["sig_q05"] = df["q_primary"] < 0.05
    df["large_effect"] = df["sig_q05"] & (df["nCG"] >= 3) & (df["diff"].abs() > 0.20)
    df["direction"] = np.where(df["diff"] > 0, "FAR_high",
                      np.where(df["diff"] < 0, "FAR_low", "none"))

    keep_cols = ["region_id", "chrom", "start", "end", "length", "nCG",
                 "control_mean", "far_mean", "diff", "direction",
                 "n_valid_control", "n_valid_far",
                 "logFC_primary", "t_primary", "p_primary", "q_primary", "df_residual_primary",
                 "logFC_sensitivity", "t_sensitivity", "p_sensitivity", "q_sensitivity",
                 "sig_q05", "large_effect",
                 "feature_class", "cpg_island_class", "nearest_gene", "distance_to_TSS",
                 "n_overlapping_genes", "overlapping_genes", "mean_cov_region"]
    ann = df[keep_cols].copy()
    ann.insert(0, "assembly", "mm39")
    ann.to_csv(os.path.join(TAB, f"DMR_annotation_{tissue}.tsv.gz"), sep="\t",
               index=False, compression="gzip", float_format="%.6g")

    sub = ann[ann["sig_q05"]].sort_values(["q_primary", "diff"])
    sub.to_csv(os.path.join(TAB, f"regions_q05_{tissue}.tsv"), sep="\t", index=False,
               float_format="%.6g")
    le = ann[ann["large_effect"]].sort_values(["q_primary", "diff"])
    le.to_csv(os.path.join(TAB, f"regions_large_effect_{tissue}.tsv"), sep="\t",
              index=False, float_format="%.6g")
    multi = ann[ann["n_overlapping_genes"] >= 2][["region_id", "chrom", "start", "end",
               "n_overlapping_genes", "overlapping_genes", "sig_q05", "q_primary", "diff"]]
    multi.to_csv(os.path.join(TAB, f"DMR_multiGene_{tissue}.tsv.gz"), sep="\t",
                 index=False, compression="gzip", float_format="%.6g")

    dss = pd.read_csv(os.path.join(TAB, f"DSS_regions_{tissue}.tsv"), sep="\t")
    dml = pd.read_csv(os.path.join(TAB, f"dss_run_qc_{tissue}.csv"))
    # DSS region vs significant region overlap
    n_dss_overlap = 0
    for _, r in dss.iterrows():
        m = (ann["chrom"] == r["chr"]) & (ann["start"] < r["end"]) & (ann["end"] > r["start"])
        n_dss_overlap += int(ann.loc[m, "sig_q05"].any()) if m.any() else 0
    summary_rows.append(dict(
        tissue=tissue, assembly="mm39",
        n_candidate_regions=len(ann),
        n_tested=int(ann["q_primary"].notna().sum()),
        n_q05_primary=int(ann["sig_q05"].sum()),
        n_q05_FARhigh=int(((ann["sig_q05"]) & (ann["diff"] > 0)).sum()),
        n_q05_FARlow=int(((ann["sig_q05"]) & (ann["diff"] < 0)).sum()),
        n_large_effect=int(ann["large_effect"].sum()),
        n_large_effect_FARhigh=int(((ann["large_effect"]) & (ann["diff"] > 0)).sum()),
        n_large_effect_FARlow=int(((ann["large_effect"]) & (ann["diff"] < 0)).sum()),
        n_q05_sensitivity=int((ann["q_sensitivity"] < 0.05).sum()),
        n_large_effect_sensitivity=int(((ann["q_sensitivity"] < 0.05) &
                                        (ann["diff"].abs() > 0.20)).sum()),
        min_q_primary=float(ann["q_primary"].min()),
        min_q_sensitivity=float(ann["q_sensitivity"].min()),
        n_dss_dml=int(dml["n_dml"].iloc[0]),
        n_dss_dml_fdr05=int(dml["n_dml_fdr05"].iloc[0]),
        n_dss_regions=int(dml["n_dss_regions"].iloc[0]),
        n_dss_regions_overlapping_q05_regions=int(n_dss_overlap),
        n_regions_multi_gene=int((ann["n_overlapping_genes"] >= 2).sum()),
        n_regions_with_gene=int((ann["n_overlapping_genes"] >= 1).sum()),
    ))
    print(pd.DataFrame(summary_rows[-1:]).T.to_string(), flush=True)

    # distribution tables
    for cls_col in ("feature_class", "cpg_island_class", "chrom", "direction"):
        a = ann[cls_col].value_counts().rename("all_regions")
        b = ann.loc[ann["sig_q05"], cls_col].value_counts().rename("q_lt_0.05")
        c = ann.loc[ann["large_effect"], cls_col].value_counts().rename("large_effect")
        d = pd.concat([a, b, c], axis=1).fillna(0).astype(int)
        d.index.name = cls_col
        d.to_csv(os.path.join(TAB, f"distribution_{cls_col}_{tissue}.csv"))
        print(f"\n[{tissue}] {cls_col}:\n", d.to_string(), flush=True)

    # delta / length distributions for figures
    dist_rows.append(dict(tissue=tissue, n=len(ann),
                          diff_mean=ann["diff"].mean(), diff_sd=ann["diff"].std(),
                          diff_q05=ann.loc[ann["sig_q05"], "diff"].mean() if ann["sig_q05"].any() else np.nan))

pd.DataFrame(summary_rows).to_csv(os.path.join(TAB, "result_summary.csv"), index=False)
print("\n=== RESULT SUMMARY ===")
print(pd.DataFrame(summary_rows).T.to_string())
