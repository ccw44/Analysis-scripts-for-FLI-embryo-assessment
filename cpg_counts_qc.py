#!/usr/bin/env python
"""
Step 01 - build per-sample M/U caches, sample QC, and the locked comparable-CpG sets.

Input : 12 x 6-column bed files (chrom, start0, end, beta, M, U) in the project root.
Output: tmp/<sample>.idx.npy, tmp/<sample>.M.npy, tmp/<sample>.U.npy
        tables/sample_QC.csv, tables/coverage_filter_summary.csv,
        tables/comparable_<tissue>.tsv.gz  (chrom,start,end + per-sample M/cov)
        tables/global_methylation.csv
        logs/01_build_and_qc.log

Locked rules (fixed before any differential result was inspected):
  comparable CpG  : coverage (M+U) >= MINCOV in at least MIN_SAMPLES of the 3 samples of
                    EACH group, per tissue.
  sensitivity set : same rule with MINCOV10.
"""
import os, sys, time, gzip, json
import numpy as np
import pandas as pd

# ===========================================================================
# Data locations - edit this block only.  Relative paths resolve against the
# working directory, so run the script from the repository root.
#   IN_CPG_DIR  per-sample CpG tables (steps 1 and 2 only), one file per sample:
#               chrom, start0, end, beta, M, U  -- tab-separated, no header
#   REF_DIR     mm39 annotation: mm39.chrom.sizes, refGene.txt.gz,
#               cpgIslandExt.txt.gz
#   OUT_DIR     output root: intermediate tables in tables/, figures in figures/
# ===========================================================================
IN_CPG_DIR = "data/cpg_tables"
REF_DIR = "data/mm39"
OUT_DIR = "results"
TAB = os.path.join(OUT_DIR, "tables")
FIG = os.path.join(OUT_DIR, "figures")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
TMP = os.path.join(OUT_DIR, "tmp")
LOG_DIR = os.path.join(OUT_DIR, "logs")
os.makedirs(TMP, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

MINCOV = 5
MIN_SAMPLES = 2
MINCOV10 = 10

SAMPLES = [
    # sample_id,        file,                  tissue,  group
    ("Control-liver-1", "Control-liver-1.bed", "liver",  "Control"),
    ("Control-liver-2", "Control-liver-2.bed", "liver",  "Control"),
    ("Control-liver-3", "Control-liver-3.bed", "liver",  "Control"),
    ("FAR-liver-1",     "FAR-liver-1.bed",     "liver",  "FAR"),
    ("FAR-liver-2",     "FAR-liver-2.bed",     "liver",  "FAR"),
    ("FAR-liver-3",     "FAR-liver-3.bed",     "liver",  "FAR"),
    ("Control-Blood-1", "Control-Blood-1.bed", "blood",  "Control"),
    ("Control-Blood-2", "Control-Blood-2.bed", "blood",  "Control"),
    ("Control-Blood-3", "Control-Blood-3.bed", "blood",  "Control"),
    ("FAR-Blood-1",     "FAR-Blood-1.bed",     "blood",  "FAR"),
    ("FAR-Blood-2",     "FAR-Blood-2.bed",     "blood",  "FAR"),
    ("FAR-Blood-3",     "FAR-Blood-3.bed",     "blood",  "FAR"),
]

logf = open(os.path.join(LOG_DIR, "cpg_qc.log"), "w")
def log(*a):
    m = " ".join(str(x) for x in a)
    print(m, flush=True); logf.write(m + "\n"); logf.flush()

def tic(): return time.time()

# ---- chromosome frame (UCSC mm39 chrom.sizes) -------------------------------
sizes = pd.read_csv(os.path.join(REF_DIR, "mm39.chrom.sizes"), sep="\t", header=None,
                    names=["chrom", "size"])
sizes = sizes.set_index("chrom")["size"].to_dict()
ORDER = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"]
OFF = {}
_cum = 0
for c in ORDER:
    OFF[c] = _cum
    _cum += sizes[c]
GENOME = _cum
OFF["chrMT"] = OFF["chrM"]          # source files use chrMT for the mitochondrion
log(f"genome frame: {len(ORDER)} sequences, cumulative length = {GENOME:,} bp")
log(f"mm39 chrM length = {sizes['chrM']}")

def load_sample(fn):
    """read one 6-column bed, return (idx int64, M int32, U int32) sorted by idx"""
    df = pd.read_csv(os.path.join(IN_CPG_DIR, fn), sep="\t", header=None, usecols=[0, 1, 2, 3, 4, 5],
                     dtype={1: "int64", 2: "int64", 3: "float64", 4: "int32", 5: "int32"},
                     engine="c")
    chrom = df[0].to_numpy()
    start = df[1].to_numpy()
    beta = df[3].to_numpy()
    M = df[4].to_numpy(); U = df[5].to_numpy()
    end = df[2].to_numpy()
    # structural checks
    assert np.all(end - start == 1), f"{fn}: non 1-bp intervals"
    cov = M.astype(np.int64) + U.astype(np.int64)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(cov > 0, M / np.where(cov > 0, cov, 1), np.nan)
    dev = np.nanmax(np.abs(ratio - beta))
    offs = np.array([OFF.get(c, -1) for c in pd.unique(chrom)])
    uniq_chrom = list(pd.unique(chrom))
    assert (offs >= 0).all(), f"{fn}: unknown chromosome {[c for c,o in zip(uniq_chrom,offs) if o<0]}"
    off_map = dict(zip(uniq_chrom, offs))
    idx = np.fromiter((off_map[c] + s for c, s in zip(chrom, start)), dtype=np.int64, count=len(chrom))
    # per-chromosome-block monotonicity (blocks themselves are in lexicographic order)
    for c in uniq_chrom:
        m = chrom == c
        assert np.all(np.diff(start[m]) > 0), f"{fn}: {c} positions not strictly increasing"
    assert idx.max() < GENOME, f"{fn}: coordinate beyond genome frame"
    # sort into the canonical genome frame (chr1..chr19, chrX, chrY, chrM) and prove uniqueness
    o = np.argsort(idx, kind="stable")
    idx, M, U = idx[o], M[o], U[o]
    d = np.diff(idx)
    assert np.all(d > 0), f"{fn}: {int((d == 0).sum())} duplicated CpG coordinates"
    meta = dict(file=fn, n_cpg=int(len(idx)), max_abs_beta_dev=float(dev),
                n_chrom=int(len(uniq_chrom)), chroms=";".join(map(str, uniq_chrom)),
                block_order=";".join(map(str, uniq_chrom)),
                genome_span_bp=int(idx.max() - idx.min()))
    return idx, M, U, meta

# ---- pass 1: per-sample cache + sample QC ----------------------------------
qc_rows, cache = [], {}
from functools import reduce
def union_all(arrs):
    return reduce(np.union1d, arrs)

import glob
for sid, fn, tissue, group in SAMPLES:
    t0 = tic()
    if os.environ.get("REUSE_CACHE") == "1" and os.path.exists(os.path.join(TMP, f"{sid}.idx.npy")):
        idx = np.load(os.path.join(TMP, f"{sid}.idx.npy"))
        M = np.load(os.path.join(TMP, f"{sid}.M.npy"))
        U = np.load(os.path.join(TMP, f"{sid}.U.npy"))
        meta = dict(file=fn, n_cpg=int(len(idx)), max_abs_beta_dev=float("nan"),
                    n_chrom=int(len(np.unique(idx // 10**9))), chroms="", block_order="",
                    genome_span_bp=int(idx.max() - idx.min()))
        reused = True
    else:
        idx, M, U, meta = load_sample(fn)
        reused = False
    if not reused:
        np.save(os.path.join(TMP, f"{sid}.idx.npy"), idx)
        np.save(os.path.join(TMP, f"{sid}.M.npy"), M)
        np.save(os.path.join(TMP, f"{sid}.U.npy"), U)
    cov = M.astype(np.int64) + U
    n = len(cov)
    row = dict(sample_id=sid, tissue=tissue, group=group, file=fn,
               n_cpg=n, n_chrom=meta["n_chrom"],
               max_abs_beta_minus_M_over_MplusU=meta["max_abs_beta_dev"],
               cov_mean=round(float(cov.mean()), 4),
               cov_median=float(np.median(cov)),
               cov_p25=float(np.percentile(cov, 25)), cov_p75=float(np.percentile(cov, 75)),
               cov_p99_9=float(np.percentile(cov, 99.9)), cov_max=int(cov.max()),
               n_cov_ge1=int(n),
               n_cov_ge5=int((cov >= 5).sum()), frac_cov_ge5=round(float((cov >= 5).mean()), 6),
               n_cov_ge10=int((cov >= 10).sum()), frac_cov_ge10=round(float((cov >= 10).mean()), 6),
               n_cov_ge20=int((cov >= 20).sum()), frac_cov_ge20=round(float((cov >= 20).mean()), 6),
               n_cov_1_2=int(((cov >= 1) & (cov <= 2)).sum()),
               frac_cov_1_2=round(float(((cov >= 1) & (cov <= 2)).mean()), 6),
               n_methylated_total=int(M.sum(dtype=np.int64)),
               n_unmethylated_total=int(U.sum(dtype=np.int64)),
               global_beta_unweighted=round(float(np.nanmean(M / cov)), 6),
               global_beta_weighted=round(float(M.sum(dtype=np.int64) / cov.sum()), 6))
    qc_rows.append(row)
    cache[sid] = dict(idx=idx, M=M, U=U, cov=cov, tissue=tissue, group=group)
    log(f"[{sid}] n={n:,} cov_mean={row['cov_mean']} ge5={row['frac_cov_ge5']:.4f} "
        f"ge10={row['frac_cov_ge10']:.4f} beta_dev={meta['max_abs_beta_dev']:.3g} "
        f"({tic()-t0:.1f}s)")
    del idx, M, U

qc = pd.DataFrame(qc_rows)
qc.to_csv(os.path.join(TAB, "sample_QC.csv"), index=False)

# ---- pass 2: per-tissue comparable sets ------------------------------------
def comparable_positions(sample_ids, mincov):
    """positions where >=MIN_SAMPLES samples of this group have cov >= mincov"""
    unions = None
    for sid in sample_ids:
        p = cache[sid]["idx"][cache[sid]["cov"] >= mincov]
        unions = p if unions is None else np.union1d(unions, p)
    cnt = np.zeros(len(unions), dtype=np.int8)
    for sid in sample_ids:
        idx = cache[sid]["idx"]
        pos = np.searchsorted(idx, unions)
        ok = (pos < len(idx)) & (idx[np.clip(pos, 0, len(idx) - 1)] == unions)
        good = np.zeros(len(unions), dtype=bool)
        good[ok] = cache[sid]["cov"][pos[ok]] >= mincov
        cnt += good.astype(np.int8)
    return unions[cnt >= MIN_SAMPLES]

def beta_matrix(positions, sample_ids):
    """per-sample M / cov / beta on `positions` (NaN where the CpG is not reported)"""
    Mmat = np.zeros((len(sample_ids), len(positions)), dtype=np.int32)
    Cmat = np.zeros((len(sample_ids), len(positions)), dtype=np.int32)
    for i, sid in enumerate(sample_ids):
        idx = cache[sid]["idx"]
        pos = np.searchsorted(idx, positions)
        pos_c = np.clip(pos, 0, len(idx) - 1)
        ok = (pos < len(idx)) & (idx[pos_c] == positions)
        Mmat[i, ok] = cache[sid]["M"][pos_c[ok]]
        Cmat[i, ok] = cache[sid]["cov"][pos_c[ok]]
    return Mmat, Cmat

cov_filter_rows, comp_sets = [], {}
for tissue in ("liver", "blood"):
    t0 = tic()
    ctrl = [s for s, f, t, g in SAMPLES if t == tissue and g == "Control"]
    far = [s for s, f, t, g in SAMPLES if t == tissue and g == "FAR"]
    tids = ctrl + far
    total_union = union_all([cache[s]["idx"] for s in tids])
    for mincov, tag in ((MINCOV, "primary_ge5"), (MINCOV10, "sensitivity_ge10")):
        pc = comparable_positions(ctrl, mincov)
        pf = comparable_positions(far, mincov)
        comp = np.intersect1d(pc, pf)
        log(f"[{tissue}/{tag}] per-group union ctrl={len(pc):,} far={len(pf):,} "
            f"comparable={len(comp):,} ({100*len(comp)/len(total_union):.2f}% of union-{tissue})")
        cov_filter_rows.append(dict(
            tissue=tissue, min_coverage=mincov, min_samples_per_group=MIN_SAMPLES,
            n_cpg_union_any_sample=len(total_union),
            n_cpg_ctrl_comparable=len(pc), n_cpg_far_comparable=len(pf),
            n_cpg_comparable=len(comp),
            frac_of_union=round(len(comp) / len(total_union), 6),
            frac_of_sample_mean_cpg=round(len(comp) / np.mean([cache[s]["idx"].size for s in tids]), 6),
            loss_from_union=len(total_union) - len(comp)))
        comp_sets[(tissue, tag)] = comp
    # write the primary comparable table for R / DSS
    comp = comp_sets[(tissue, "primary_ge5")]
    Mmat, Cmat = beta_matrix(comp, tids)
    chroms = np.array(ORDER)
    # decode idx -> chrom/start using offsets
    offs = np.array([OFF[c] for c in ORDER[:-1]] + [OFF["chrM"]])
    ci = np.searchsorted(offs, comp, side="right") - 1
    starts = comp - offs[ci]
    out = pd.DataFrame({"chrom": chroms[ci], "start": starts, "end": starts + 1})
    for i, sid in enumerate(tids):
        out[sid + "_M"] = Mmat[i]
        out[sid + "_cov"] = Cmat[i]
    fn = os.path.join(TAB, f"comparable_{tissue}.tsv.gz")
    out.to_csv(fn, sep="\t", index=False, compression="gzip")
    log(f"[{tissue}] wrote {fn} ({out.shape[0]:,} rows x {out.shape[1]} cols, {tic()-t0:.1f}s)")
    # effective CpG spacing on the comparable set
    sp = np.diff(comp)
    sp = sp[sp > 0]
    log(f"[{tissue}] comparable CpG spacing: median={np.median(sp):.0f}bp "
        f"frac<100bp={np.mean(sp<100):.3f}")
    del Mmat, Cmat, out

pd.DataFrame(cov_filter_rows).to_csv(os.path.join(TAB, "coverage_filter_summary.csv"), index=False)

# ---- global methylation table ---------------------------------------------
gm = []
for tissue in ("liver", "blood"):
    ctrl = [s for s, f, t, g in SAMPLES if t == tissue and g == "Control"]
    far = [s for s, f, t, g in SAMPLES if t == tissue and g == "FAR"]
    comp = comp_sets[(tissue, "primary_ge5")]
    Mmat, Cmat = beta_matrix(comp, ctrl + far)
    with np.errstate(invalid="ignore", divide="ignore"):
        b = np.where(Cmat > 0, Mmat / np.where(Cmat > 0, Cmat, 1), np.nan)
    for i, sid in enumerate(ctrl + far):
        v = b[i]
        gm.append(dict(sample_id=sid, tissue=tissue,
                       group="Control" if sid in ctrl else "FAR",
                       n_cpg_used=int(np.isfinite(v).sum()),
                       mean_beta_unweighted=float(np.nanmean(v)),
                       sd_beta_unweighted=float(np.nanstd(v)),
                       mean_beta_weighted=float(Mmat[i].sum() / Cmat[i].sum()),
                       mean_cov=float(Cmat[i].mean())))
    log(f"[{tissue}] global beta: " + ", ".join(
        f"{r['sample_id']}={r['mean_beta_unweighted']:.4f}" for r in gm if r["tissue"] == tissue))
    del Mmat, Cmat, b
pd.DataFrame(gm).to_csv(os.path.join(TAB, "global_methylation.csv"), index=False)

# ---- sample correlation / PCA on all-sample-complete CpGs ------------------
from numpy.linalg import svd
rel_rows = []
for tissue in ("liver", "blood"):
    ctrl = [s for s, f, t, g in SAMPLES if t == tissue and g == "Control"]
    far = [s for s, f, t, g in SAMPLES if t == tissue and g == "FAR"]
    tids = ctrl + far
    comp = comp_sets[(tissue, "primary_ge5")]
    Mmat, Cmat = beta_matrix(comp, tids)
    complete = np.all(Cmat >= MINCOV, axis=0)
    log(f"[{tissue}] CpGs with cov>=5 in ALL {len(tids)} samples: {complete.sum():,} "
        f"({100*complete.mean():.2f}% of comparable set)")
    sub = comp[complete]
    Ms, Cs = Mmat[:, complete], Cmat[:, complete]
    b = Ms / Cs
    # pairwise complete correlation on the comparable set (primary), and on the complete set
    with np.errstate(invalid="ignore"):
        B = np.where(Cmat > 0, Mmat / np.where(Cmat > 0, Cmat, 1), np.nan)
    R = np.ma.corrcoef(np.ma.masked_invalid(B)).filled(np.nan)
    Rc = np.corrcoef(b)
    for i, si in enumerate(tids):
        for j, sj in enumerate(tids):
            rel_rows.append(dict(tissue=tissue, sample_i=si, sample_j=sj,
                                 pearson_comparableSet=round(float(R[i, j]), 5),
                                 pearson_completeSet=round(float(Rc[i, j]), 5)))
    pd.DataFrame(R, index=tids, columns=tids).to_csv(
        os.path.join(TAB, f"sample_correlation_comparable_{tissue}.csv"),
        float_format="%.5f")
    pd.DataFrame(Rc, index=tids, columns=tids).to_csv(
        os.path.join(TAB, f"sample_correlation_complete_{tissue}.csv"),
        float_format="%.5f")
    # PCA on the complete set, centred per CpG
    X = b - b.mean(axis=1, keepdims=True)
    # PCA across samples: samples x CpGs -> use correlation-free svd on the 6 x n matrix
    U_, s_, Vt_ = svd(X, full_matrices=False)
    ev = s_ ** 2
    pcs = U_ * s_
    pca = pd.DataFrame(pcs, index=tids, columns=[f"PC{k+1}" for k in range(pcs.shape[1])])
    pca.insert(0, "group", ["Control" if t in ctrl else "FAR" for t in tids])
    pca.to_csv(os.path.join(TAB, f"sample_PCA_{tissue}.csv"), float_format="%.5f")
    pd.Series(ev / ev.sum(), index=[f"PC{k+1}" for k in range(len(ev))]).to_csv(
        os.path.join(TAB, f"sample_PCA_varexpl_{tissue}.csv"), header=["variance_explained"])
    log(f"[{tissue}] PCA var explained: " + ", ".join(f"PC{i+1}={v:.3f}" for i, v in enumerate(ev / ev.sum())))
    del Mmat, Cmat, B, b, X

pd.DataFrame(rel_rows).to_csv(os.path.join(TAB, "sample_correlation_pairs.csv"), index=False)
log("DONE 01")
logf.close()
