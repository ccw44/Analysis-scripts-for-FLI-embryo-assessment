#!/usr/bin/env python
"""Candidate regions, per-sample region methylation and the effect size.

Geometry is fixed from coverage only: merge comparable CpGs with a gap <= GAP_MAX bp and
require >= MIN_CG CpGs. A sample is valid for a region when it covers >= 50% of the
region CpGs and at least two CpGs; the group mean is the equal-weight mean over its valid
samples, and the effect size is FLI 2 h minus control.
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


GAP_MAX = 100
MIN_CG = 3
MIN_FRAC_COVERED = 0.50
MIN_N_COVERED = 2
MIN_VALID_PER_GROUP = 2

ORDER = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"]
CIDX = {c: i for i, c in enumerate(ORDER)}

for tissue, ctrl, far in (
    ("liver", ["Control-liver-1", "Control-liver-2", "Control-liver-3"],
              ["FAR-liver-1", "FAR-liver-2", "FAR-liver-3"]),
    ("blood", ["Control-Blood-1", "Control-Blood-2", "Control-Blood-3"],
              ["FAR-Blood-1", "FAR-Blood-2", "FAR-Blood-3"]),
):
    print(f"\n=== {tissue} ===", flush=True)
    dt = pd.read_csv(os.path.join(TAB, f"comparable_{tissue}.tsv.gz"), sep="\t",
                     dtype={"chrom": "string", "start": np.int64, "end": np.int64})
    sids = far + ctrl
    chrom = dt["chrom"].to_numpy()
    start = dt["start"].to_numpy()
    cc = np.array([CIDX.get(c, -1) for c in chrom])
    assert (cc >= 0).all(), "unknown chromosome"

    # ---- region breaks: new chromosome, or gap > GAP_MAX --------------------
    brk = np.empty(len(start), dtype=bool)
    brk[0] = True
    brk[1:] = (cc[1:] != cc[:-1]) | ((start[1:] - start[:-1]) > GAP_MAX)
    region_id = np.cumsum(brk) - 1
    n_reg_raw = region_id[-1] + 1

    # start/end/nCG per region
    first = np.flatnonzero(brk)
    reg_chrom = cc[first]
    reg_start = start[first]
    last = np.append(first[1:] - 1, len(start) - 1)
    reg_end = dt["end"].to_numpy()[last]
    reg_ncg = last - first + 1

    keep = reg_ncg >= MIN_CG
    print(f"raw runs {n_reg_raw:,}  ->  regions with >={MIN_CG} CpG: {keep.sum():,} "
          f"({100*keep.mean():.1f}%)", flush=True)

    # ---- per-sample per-region aggregates ----------------------------------
    keep_cpg = keep[region_id]                 # CpG-level mask (drop short regions)
    cov = dt[[s + "_cov" for s in sids]].to_numpy(dtype=np.int32)[keep_cpg]
    M = dt[[s + "_M" for s in sids]].to_numpy(dtype=np.int32)[keep_cpg]
    sidx = region_id[keep_cpg]                 # CpG -> kept region id (contiguous)
    new_id = np.cumsum(np.r_[True, np.diff(sidx) > 0]) - 1
    n_reg = new_id[-1] + 1
    rows_first = np.r_[True, np.diff(new_id) > 0]
    f2 = np.flatnonzero(rows_first)
    l2 = np.r_[f2[1:] - 1, len(sidx) - 1]
    ncg = (l2 - f2 + 1)
    reg_chrom = reg_chrom[keep]; reg_start = reg_start[keep]; reg_end = reg_end[keep]

    beta = np.full((n_reg, len(sids)), np.nan, dtype=np.float32)
    ncov = np.zeros((n_reg, len(sids)), dtype=np.int32)
    meancov = np.zeros((n_reg, len(sids)), dtype=np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        b_cpg = np.where(cov > 0, M / np.where(cov > 0, cov, 1), np.nan)
    for j, s in enumerate(sids):
        cj = cov[:, j]; bj = b_cpg[:, j]
        covd = (cj > 0).astype(np.float64)
        val = np.where(cj > 0, bj, 0.0).astype(np.float64)
        sc = np.add.reduceat(covd, f2)
        sv = np.add.reduceat(val, f2)
        ssum_cov = np.add.reduceat(cj.astype(np.float64), f2)
        ncov[:, j] = sc.astype(np.int32)
        with np.errstate(invalid="ignore", divide="ignore"):
            beta[:, j] = np.where(sc > 0, sv / np.where(sc > 0, sc, 1), np.nan)
            meancov[:, j] = np.where(sc > 0, ssum_cov / np.where(sc > 0, sc, 1), np.nan)

    # ---- eligibility + contract effect size --------------------------------
    need = np.maximum(MIN_N_COVERED, np.ceil(MIN_FRAC_COVERED * ncg)).astype(np.int32)
    valid = ncov >= need[:, None]
    nv_c = valid[:, [sids.index(c) for c in ctrl]].sum(axis=1)
    nv_f = valid[:, [sids.index(f) for f in far]].sum(axis=1)

    def group_mean(cols):
        idx = [sids.index(c) for c in cols]
        v = valid[:, idx]; b = beta[:, idx]
        with np.errstate(invalid="ignore"):
            s = np.where(v, b, 0).sum(axis=1, dtype=np.float64)
            n = v.sum(axis=1)
            m = np.where(n > 0, s / np.where(n > 0, n, 1), np.nan)
        return m, n

    ctrl_mean, n_c = group_mean(ctrl)
    far_mean, n_f = group_mean(far)
    evaluable = (n_c >= MIN_VALID_PER_GROUP) & (n_f >= MIN_VALID_PER_GROUP)
    diff = np.where(evaluable, far_mean - ctrl_mean, np.nan)
    testable = (n_c + n_f) >= 3        # >= 3 valid samples total for a residual df >= 1

    out = pd.DataFrame({
        "region_id": np.arange(n_reg),
        "chrom": [ORDER[c] for c in reg_chrom],
        "start": reg_start, "end": reg_end,
        "length": reg_end - reg_start,
        "nCG": ncg,
        "n_valid_control": n_c, "n_valid_far": n_f,
        "control_mean": ctrl_mean, "far_mean": far_mean, "diff": diff,
        "effect_evaluable": evaluable, "testable": testable,
    })
    for j, s in enumerate(sids):
        out[s + "_beta"] = beta[:, j]
        out[s + "_nCpG_covered"] = ncov[:, j]
        out[s + "_mean_cov"] = meancov[:, j]

    fn = os.path.join(TAB, f"regions_all_{tissue}.tsv.gz")
    out.to_csv(fn, sep="\t", index=False, compression="gzip", float_format="%.6f")
    mm = out[["region_id"]].copy()
    for s in sids:
        mm[s] = np.where(valid[:, sids.index(s)], beta[:, sids.index(s)], np.nan)
    mm.to_csv(os.path.join(TAB, f"region_mean_matrix_{tissue}.tsv.gz"), sep="\t",
              index=False, compression="gzip", float_format="%.6f")

    geo = pd.DataFrame([dict(
        tissue=tissue, gap_max_bp=GAP_MAX, min_cpg=MIN_CG,
        n_comparable_cpg=int(len(dt)), n_raw_runs=int(n_reg_raw),
        n_regions=int(n_reg), n_regions_ge3cpg=int(keep.sum()),
        n_regions_evaluable=int(evaluable.sum()), n_regions_testable=int(testable.sum()),
        frac_evaluable=round(float(evaluable.mean()), 4),
        median_length=int(np.median(reg_end - reg_start)),
        p90_length=int(np.percentile(reg_end - reg_start, 90)),
        p99_length=int(np.percentile(reg_end - reg_start, 99)),
        max_length=int((reg_end - reg_start).max()),
        median_nCG=int(np.median(ncg)), max_nCG=int(ncg.max()),
        frac_regions_len_gt_2000=round(float(np.mean((reg_end - reg_start) > 2000)), 4),
        frac_regions_len_gt_5000=round(float(np.mean((reg_end - reg_start) > 5000)), 4))])
    geo.to_csv(os.path.join(TAB, f"region_geometry_{tissue}.csv"), index=False)
    print(geo.T.to_string(), flush=True)
    print("effect size distribution (evaluable):\n",
          out.loc[evaluable, "diff"].describe().to_string(), flush=True)
    print(f"wrote {fn}", flush=True)
    del dt, cov, M, b_cpg, beta, ncov, meancov, out, mm
