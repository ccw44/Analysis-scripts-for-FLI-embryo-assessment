#!/usr/bin/env python
"""
Step 02 - independent re-verification of the column semantics of all 12 source files.

For every file, streamed in 2M-row chunks (never the whole file in memory):
  * beta (col 4) == M/(M+U) exactly, max absolute deviation
  * M (col 5) and U (col 6) are non-negative integers
  * M+U > 0 everywhere (no undefined ratio)
  * end - start == 1 (single-base CpG resolution)
  * coordinates inside the mm39 chromosome frame
Writes tables/columns_check.csv and fills the beta-deviation column of sample_QC.csv.
"""
import os
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

GROUPS = [("Control", "liver", 1), ("Control", "liver", 2), ("Control", "liver", 3),
          ("FAR", "liver", 1), ("FAR", "liver", 2), ("FAR", "liver", 3),
          ("Control", "Blood", 1), ("Control", "Blood", 2), ("Control", "Blood", 3),
          ("FAR", "Blood", 1), ("FAR", "Blood", 2), ("FAR", "Blood", 3)]
sizes = pd.read_csv(os.path.join(REF_DIR, "mm39.chrom.sizes"), sep="\t", header=None,
                    names=["chrom", "size"]).set_index("chrom")["size"].to_dict()
sizes["chrMT"] = sizes["chrM"]

rows = []
for g, tissue_label, k in GROUPS:
    sid = f"{g}-{tissue_label}-{k}"
    fn = os.path.join(IN_CPG_DIR, sid + ".bed")
    max_dev = 0.0
    n = 0
    nneg = 0
    nzero = 0
    nwidth = 0
    oob = 0
    for chunk in pd.read_csv(fn, sep="\t", header=None, usecols=[0, 1, 2, 3, 4, 5],
                             dtype={1: "int64", 2: "int64", 3: "float64", 4: "int64", 5: "int64"},
                             chunksize=2_000_000, engine="c"):
        chrom, start, end, beta, M, U = (chunk[c].to_numpy() for c in range(6))
        cov = M + U
        nzero += int((cov == 0).sum())
        nneg += int(((M < 0) | (U < 0)).sum())
        nwidth += int((end - start != 1).sum())
        if len(chrom):
            lim = np.array([sizes.get(c, -1) for c in chrom])
            oob += int(((start < 0) | (start >= lim)).sum())
        ok = cov > 0
        if ok.any():
            dev = np.abs(beta[ok] - M[ok] / cov[ok])
            max_dev = max(max_dev, float(dev.max()))
        n += len(chunk)
    rows.append(dict(sample_id=sid, n_rows=n,
                     max_abs_beta_minus_M_over_MplusU=max_dev,
                     n_zero_coverage=nzero, n_negative_counts=nneg,
                     n_nonunit_interval=nwidth, n_out_of_frame=oob))
    print(f"{sid}: n={n:,} max_dev={max_dev:.3g} zero_cov={nzero} neg={nneg} "
          f"nonunit={nwidth} oob={oob}", flush=True)

chk = pd.DataFrame(rows)
chk.to_csv(os.path.join(TAB, "columns_check.csv"), index=False)

qc = pd.read_csv(os.path.join(TAB, "sample_QC.csv"))
qc = qc.drop(columns=["max_abs_beta_minus_M_over_MplusU"]).merge(
    chk[["sample_id", "max_abs_beta_minus_M_over_MplusU"]], on="sample_id", how="left")
qc.to_csv(os.path.join(TAB, "sample_QC.csv"), index=False)
print("\nALL CHECKS: max deviation across all files =",
      chk["max_abs_beta_minus_M_over_MplusU"].max())
print("total zero-coverage rows:", chk["n_zero_coverage"].sum(),
      "| negative counts:", chk["n_negative_counts"].sum(),
      "| non-unit intervals:", chk["n_nonunit_interval"].sum(),
      "| out-of-frame:", chk["n_out_of_frame"].sum())
