#!/usr/bin/env python
"""All methylation figures: a vector PDF per panel plus a CSV of the plotted values.

Data locations are set in the block below; run from the repository root.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.ticker import FuncFormatter, MaxNLocator

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


# ---- style: OmicVerse if available, else a matching matplotlib fallback -----
avail = {f.name for f in fm.fontManager.ttflist}
sans = next((f for f in ("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans")
             if f in avail), "DejaVu Sans")
try:
    import omicverse as ov
    fp = fm.findfont(fm.FontProperties(family=sans), fallback_to_default=False)
    ov.style(font_path=fp)
    RENDERER = "omicverse ov.style + matplotlib"
except Exception as e:
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [sans, "DejaVu Sans"],
                         "axes.linewidth": 0.8, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42,
                         "ps.fonttype": 42, "svg.fonttype": "none"})
    RENDERER = f"matplotlib fallback ({type(e).__name__}: {e})"
print("renderer:", RENDERER)

import plot_style as _style      # single font family, no per-glyph fallback
FAMILY = _style.apply()
print("font family:", FAMILY)

COL_SAMPLE = {"Control": "#0072B2", "FAR": "#D55E00"}
CB5 = ["#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2"]
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 300})


def save(fig, name, df=None):
    # PDF only: vector, editable text (Type 42). No PNG is written any more.
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    p_pdf = os.path.join(FIG, name + ".pdf")
    fig.savefig(p_pdf, bbox_inches="tight")
    plt.close(fig)
    if df is not None:
        df.to_csv(os.path.join(FIG, name + ".csv"), index=False)
    print("wrote", p_pdf)


# ===========================================================================
# sample QC, sample structure, region effects
# ===========================================================================
qc = pd.read_csv(os.path.join(TAB, "sample_QC.csv"))
qc["label"] = qc["sample_id"].str.replace("-Blood-", " blood ").str.replace("-liver-", " liver ")
gm = pd.read_csv(os.path.join(TAB, "global_methylation.csv"))
gm["label"] = gm["sample_id"].str.replace("-Blood-", " blood ").str.replace("-liver-", " liver ")

# ---- sample QC ----
q = qc.copy()
q["b_1_2"] = q["n_cov_1_2"]
q["b_3_4"] = q["n_cpg"] - q["n_cov_1_2"] - q["n_cov_ge5"]
q["b_5_9"] = q["n_cov_ge5"] - q["n_cov_ge10"]
q["b_10_19"] = q["n_cov_ge10"] - q["n_cov_ge20"]
q["b_20p"] = q["n_cov_ge20"]
tot = q["n_cpg"]
for b in ("b_1_2", "b_3_4", "b_5_9", "b_10_19", "b_20p"):
    q[b + "_pct"] = 100 * q[b] / tot

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.8), gridspec_kw={"width_ratios": [1.45, 1]})
x = np.arange(len(q)); bottom = np.zeros(len(q))
for b, lab_, col in (("b_1_2", "1–2×", "#D55E00"), ("b_3_4", "3–4×", "#E69F00"),
                     ("b_5_9", "5–9×", "#56B4E9"), ("b_10_19", "10–19×", "#009E73"),
                     ("b_20p", "≥20×", "#0072B2")):
    ax1.bar(x, q[b + "_pct"], bottom=bottom, color=col, width=0.72, label=lab_,
            edgecolor="white", linewidth=0.4)
    bottom += q[b + "_pct"].values
ax1.set_xticks(x); ax1.set_xticklabels(q["label"], rotation=45, ha="right", fontsize=8)
ax1.set_ylabel("CpG (%)"); ax1.set_ylim(0, 124); ax1.set_yticks([0, 20, 40, 60, 80, 100])
ax1.set_title("A   CpG coverage distribution", loc="left", fontsize=11, fontweight="bold")
ax1.legend(title="reads per CpG", fontsize=8, title_fontsize=8, frameon=False,
           loc="upper center", ncol=5, columnspacing=1.2, handlelength=1.0, bbox_to_anchor=(0.5, 1.01))
ax1.grid(axis="y", lw=0.5, alpha=0.35); ax1.set_axisbelow(True)

CENTERS = {"liver": {"Control": 0.0, "FAR": 1.0}, "blood": {"Control": 3.0, "FAR": 4.0}}
scat_rows = []
for t in ("liver", "blood"):
    for g in ("Control", "FAR"):
        s_ = gm[(gm.tissue == t) & (gm.group == g)].reset_index(drop=True)
        xc = CENTERS[t][g]
        xs = xc + np.linspace(-0.16, 0.16, len(s_))
        ax2.scatter(xs, s_["mean_beta_unweighted"], s=46, color=COL_SAMPLE[g], zorder=3,
                    edgecolor="white", linewidth=0.7, label=g if (t == "liver") else None)
        m, sd = s_["mean_beta_unweighted"].mean(), s_["mean_beta_unweighted"].std(ddof=1)
        ax2.errorbar([xc], [m], yerr=[sd], fmt="_", color=COL_SAMPLE[g], capsize=6, lw=1.8, ms=16, zorder=2)
        scat_rows.append(pd.DataFrame({"tissue": t, "group": g, "x_center": xc,
                                       "sample_id": s_["sample_id"], "x": xs,
                                       "mean_beta": s_["mean_beta_unweighted"],
                                       "group_mean": m, "group_sd": sd}))
ax2.set_xticks([0, 1, 3, 4])
ax2.set_xticklabels(["Ctrl\nliver", "FAR\nliver", "Ctrl\nblood", "FAR\nblood"], fontsize=8.5)
ax2.set_xlim(-0.6, 4.6)
ax2.set_ylabel("Mean CpG methylation (0–1)")
ax2.set_title("B   Global methylation per sample", loc="left", fontsize=11, fontweight="bold")
ax2.legend(fontsize=8, frameon=False, loc="lower right", ncol=2)
ax2.grid(axis="y", lw=0.5, alpha=0.35); ax2.set_axisbelow(True)
fig.tight_layout()
save(fig, "fig01_sample_qc", pd.concat(scat_rows, ignore_index=True))

# ---- sample structure ----
cl = pd.read_csv(os.path.join(TAB, "sample_correlation_complete_liver.csv"), index_col=0)
cb = pd.read_csv(os.path.join(TAB, "sample_correlation_complete_blood.csv"), index_col=0)
pl = pd.read_csv(os.path.join(TAB, "sample_PCA_liver.csv"), index_col=0)
pb = pd.read_csv(os.path.join(TAB, "sample_PCA_blood.csv"), index_col=0)
vl = pd.read_csv(os.path.join(TAB, "sample_PCA_varexpl_liver.csv"), index_col=0).iloc[:, 0]
vb = pd.read_csv(os.path.join(TAB, "sample_PCA_varexpl_blood.csv"), index_col=0).iloc[:, 0]

fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.8))
pca_rows = []
for k, (mat, t) in enumerate(((cl, "liver"), (cb, "blood"))):
    ax = axes[0, k]
    vmin = min(0.93, float(np.min(mat.values)))
    # pcolormesh (not imshow) so the heatmap stays vector in the PDF
    im = ax.pcolormesh(np.arange(len(mat) + 1), np.arange(len(mat) + 1), mat.values,
                       cmap="viridis", vmin=vmin, vmax=1.0, shading="flat")
    ax.set_xlim(0, len(mat)); ax.set_ylim(0, len(mat))
    ax.invert_yaxis()
    ax.set_xticks(np.arange(len(mat)) + 0.5); ax.set_yticks(np.arange(len(mat)) + 0.5)
    lab = [c.replace("-", "\n", 1) for c in mat.index]
    ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=6.8)
    ax.set_yticklabels(lab, fontsize=6.8)
    ax.tick_params(length=0)
    for i in range(len(mat)):
        for j in range(len(mat)):
            ax.text(j + 0.5, i + 0.5, f"{mat.values[i, j]:.3f}", ha="center", va="center", fontsize=5.6,
                    color="white" if mat.values[i, j] < vmin + 0.6 * (1 - vmin) else "black")
    ax.set_title(f"{'AB'[k]}   Pearson r, complete-case CpG set ({t})", loc="left",
                 fontsize=10.5, fontweight="bold")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)

for k, (pr, v, t) in enumerate(((pl, vl, "liver"), (pb, vb, "blood"))):
    ax = axes[1, k]
    for g in ("Control", "FAR"):
        s_ = pr[pr["group"] == g]
        ax.scatter(s_["PC1"], s_["PC2"], s=78, color=COL_SAMPLE[g], alpha=0.9,
                   edgecolor="white", linewidth=0.9, label=g)
    ax.axhline(0, color="0.85", lw=0.8, zorder=0)
    ax.set_xlabel(f"PC1 ({100*v.iloc[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100*v.iloc[1]:.1f}%)")
    ax.set_title(f"{'CD'[k]}   PCA ({t}, {len(pr)} samples)", loc="left",
                 fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=8, frameon=False, loc="best")
    ax.grid(lw=0.5, alpha=0.3); ax.set_axisbelow(True)
    pca_rows.append(pr.assign(tissue=t).rename_axis("sample_id").reset_index()[
        ["sample_id", "tissue", "group", "PC1", "PC2"]])
fig.tight_layout()
save(fig, "fig02_sample_structure", pd.concat(pca_rows, ignore_index=True))

# ---- region-level effects ----
fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.2))
delta_rows = []
for j, t in enumerate(("liver", "blood")):
    d = pd.read_csv(os.path.join(TAB, f"DMR_annotation_{t}.tsv.gz"), sep="\t",
                    usecols=["region_id", "chrom", "start", "end", "diff", "q_primary",
                             "q_sensitivity", "feature_class", "cpg_island_class", "nCG"])
    d = d[d["diff"].notna()]
    ax = axes[0, j]
    ax.hist(d["diff"], bins=160, color="#56B4E9", edgecolor="none")
    ax.axvline(0, color="0.3", lw=0.8)
    ax.axvline(0.20, color="#D55E00", lw=0.8, ls="--")
    ax.axvline(-0.20, color="#D55E00", lw=0.8, ls="--")
    ax.set_xlim(-0.55, 0.55)
    ax.margins(y=0.10)
    ax.set_xlabel("Region Δ methylation (FAR − Control)")
    ax.set_ylabel("candidate regions")
    ax.set_title(f"{'AB'[j]}   {t}: region effect sizes (n={len(d):,})", loc="left",
                 fontsize=10.5, fontweight="bold")
    ax.grid(lw=0.5, alpha=0.3); ax.set_axisbelow(True)
    q = d["q_primary"].clip(lower=1e-12)
    ax = axes[1, j]
    sig = d["q_primary"] < 0.05
    ax.scatter(d.loc[~sig, "diff"], -np.log10(q[~sig]), s=3, color="0.72", alpha=0.35,
               edgecolor="none", label="q ≥ 0.05")
    if sig.any():
        ax.scatter(d.loc[sig, "diff"], -np.log10(q[sig]), s=46, color="#D55E00",
                   edgecolor="black", linewidth=0.6, label="q < 0.05", zorder=4)
    else:
        ax.text(0.02, 0.93, "no region with q < 0.05", transform=ax.transAxes,
                fontsize=8, color="#D55E00", va="top")
    ax.axhline(-np.log10(0.05), color="#0072B2", lw=0.9, ls="--")
    ax.axvline(0.20, color="0.4", lw=0.8, ls=":")
    ax.axvline(-0.20, color="0.4", lw=0.8, ls=":")
    ax.set_xlabel("Region Δ methylation (FAR − Control)")
    ax.set_ylabel("−log10 region q (BH)")
    ax.set_title(f"{'CD'[j]}   {t}: volcano", loc="left", fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax.grid(lw=0.5, alpha=0.3); ax.set_axisbelow(True)
    delta_rows.append(d.assign(tissue=t))
fig.tight_layout()
save(fig, "fig03_region_effects", pd.concat(delta_rows, ignore_index=True))

# ---- composition of the testable set ----
fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.4))
comp_rows = []
for k, (cls, title, order) in enumerate((
        ("feature_class", "A   Genomic feature class", ["promoter", "exon", "intron", "intergenic"]),
        ("cpg_island_class", "B   CpG-island context", ["island", "shore", "shelf", "open_sea"]))):
    ax = axes[k]
    w = 0.38
    for i, t in enumerate(("liver", "blood")):
        d = pd.read_csv(os.path.join(TAB, f"DMR_annotation_{t}.tsv.gz"), sep="\t", usecols=[cls])
        v = d[cls].value_counts(normalize=True).reindex(order).fillna(0) * 100
        ax.bar(np.arange(len(order)) + (i - 0.5) * w, v.values, width=w,
               color=COL_SAMPLE["Control"] if t == "liver" else COL_SAMPLE["FAR"], label=t)
        comp_rows.append(pd.DataFrame({"class_scheme": cls, "class": order,
                                       "tissue": t, "pct_of_candidate_regions": v.values}))
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=20, ha="right")
    ax.set_ylabel("% of candidate regions")
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(axis="y", lw=0.5, alpha=0.3); ax.set_axisbelow(True)
fig.tight_layout()
save(fig, "fig04_genomic_distribution", pd.concat(comp_rows, ignore_index=True))

# ---- the one significant region ----
REG = dict(chrom="chr13", start=70881700, end=70882100)
sids_c = ["Control-Blood-1", "Control-Blood-2", "Control-Blood-3"]
sids_f = ["FAR-Blood-1", "FAR-Blood-2", "FAR-Blood-3"]
dc = pd.read_csv(os.path.join(TAB, "comparable_blood.tsv.gz"), sep="\t")
sub = dc[(dc.chrom == REG["chrom"]) & (dc.start >= REG["start"]) & (dc.start < REG["end"])].copy()
cgi = pd.read_csv(os.path.join(REF_DIR, "cpgIslandExt.txt.gz"), sep="\t", header=None,
                  names=["bin", "chrom", "s", "e", "name", "length", "cpgNum", "gcNum",
                         "perCpg", "perGc", "obsExp"])
cgi = cgi[(cgi.chrom == REG["chrom"]) & (cgi.e < REG["start"]) & (cgi.e > REG["start"] - 20000)]

fig, axes = plt.subplots(2, 1, figsize=(11.0, 5.6), sharex=True,
                         gridspec_kw={"height_ratios": [3, 1]})
ax = axes[0]
rows = []
for sids, g in ((sids_c, "Control"), (sids_f, "FAR")):
    for s in sids:
        b = (sub[s + "_M"] / sub[s + "_cov"]).to_numpy()
        ax.plot(sub["start"], b, marker="o", ms=4.6, lw=1.1, color=COL_SAMPLE[g], alpha=0.9,
                label=s.replace("-Blood-", " blood "))
        rows.append(pd.DataFrame({"chrom": sub["chrom"], "start": sub["start"], "sample": s,
                                  "group": g, "beta": b, "cov": sub[s + "_cov"]}))
ax.axvspan(70881878, 70881916, color="0.55", alpha=0.22, zorder=0)
ax.set_ylabel("CpG methylation (0–1)"); ax.set_ylim(-0.06, 1.10)
ax.tick_params(axis="x", labelbottom=True)
ax.set_title("A   chr13:70,881,700–70,882,100  (DMR 70,881,878–70,881,916; q=0.0041)",
             loc="left", fontsize=10.5, fontweight="bold")
ax.grid(lw=0.4, alpha=0.3); ax.set_axisbelow(True)
ax = axes[1]
for sids, g in ((sids_c, "Control"), (sids_f, "FAR")):
    for s in sids:
        ax.plot(sub["start"], sub[s + "_cov"], marker="o", ms=3.4, lw=0.9, color=COL_SAMPLE[g], alpha=0.75)
ax.axvspan(70881878, 70881916, color="0.55", alpha=0.22, zorder=0)
ax.set_ylim(0, 11); ax.set_yticks([0, 5, 10])
ax.set_ylabel("coverage"); ax.set_xlabel(f"{REG['chrom']} position (bp, mm39)")
ax.set_title("B   per-sample coverage", loc="left", fontsize=10.5, fontweight="bold")
ax.grid(lw=0.4, alpha=0.3); ax.set_axisbelow(True)
fig.tight_layout()
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.02),
           handlelength=1.5, columnspacing=1.1, fontsize=7.5)
save(fig, "fig07_representative_region", pd.concat(rows, ignore_index=True))

print("\nrenderer used for all figures:", RENDERER)


# ===========================================================================
# region counts and chromosome distribution
# ===========================================================================
COL_TISSUE = {"liver": "#0072B2", "blood": "#D55E00"}
DIRCOL = {"FAR_high": "#D55E00", "FAR_low": "#0072B2", "none": "#9E9E9E"}

res = pd.read_csv(os.path.join(TAB, "result_summary.csv")).set_index("tissue")

# ---- counts + direction ----
dirs = {t: pd.read_csv(os.path.join(TAB, f"distribution_direction_{t}.csv")).set_index("direction")
        for t in ("liver", "blood")}
ORDER = ["FAR_high", "FAR_low", "none"]

fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6),
                         gridspec_kw={"width_ratios": [1.0, 1.25]})

# --- A: candidate regions by direction
ax = axes[0]
w = 0.36
rows_a = []
for i, t in enumerate(("liver", "blood")):
    vals = [int(dirs[t].loc[d, "all_regions"]) for d in ORDER]
    tot = sum(vals)
    ax.bar(np.arange(len(ORDER)) + (i - 0.5) * w, vals, width=w, color=COL_TISSUE[t], label=t)
    for j, v in enumerate(vals):
        ax.text(j + (i - 0.5) * w, v * 1.015, f"{v:,}\n({100*v/tot:.1f}%)", ha="center", va="bottom",
                fontsize=6.8)
    for d, v in zip(ORDER, vals):
        rows_a.append(dict(panel="A_candidate_regions_by_direction", tissue=t, direction=d,
                           count=v, pct_of_tissue=round(100 * v / tot, 4)))
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(["FAR hyper-\nmethylated", "FAR hypo-\nmethylated", "no call\n(|Δ|<0.05)"], fontsize=8.5)
ax.set_ylabel("candidate regions")
ax.set_ylim(0, 82000)
ax.set_yticks([0, 20000, 40000, 60000, 80000])
ax.set_yticklabels(["0", "20,000", "40,000", "60,000", "80,000"])
ax.set_title("A   Candidate regions by direction", loc="left", fontsize=11, fontweight="bold")
ax.legend(fontsize=8.5, frameon=False, loc="upper right")
ax.grid(axis="y", lw=0.5, alpha=0.35); ax.set_axisbelow(True)

# --- B: region-level significant counts only (no per-CpG / site-level numbers) ---
ax = axes[1]
LADDER = [
    ("DSS-called regions  (differentially methylated)", "n_dss_regions",     "region"),
    ("region q<0.05  (primary model)",      "n_q05_primary",          "region"),
    ("region q<0.05  (sensitivity model)",  "n_q05_sensitivity",      "region"),
    ("q<0.05 & |Δ|>0.20",                   "n_large_effect",         "region"),
]
ypos = np.arange(len(LADDER))[::-1]
h = 0.36
for i, t in enumerate(("liver", "blood")):
    vals = [max(res.loc[t, k], 0) for _, k, _ in LADDER]
    yy = ypos + (0.5 - i) * h
    ax.barh(yy, vals, height=h, color=COL_TISSUE[t], label=t)
    for y_, v in zip(yy, vals):
        # zeros are annotated at the axis origin in the tissue colour so the
        # label is attributable to its series even though the bar has no length
        ax.text(v + 0.34, y_, f"{int(v):,}", va="center", ha="left", fontsize=8,
                color=COL_TISSUE[t] if v == 0 else "black")
ax.set_xlim(0, 31)
ax.set_xlabel("number of regions")
ax.set_yticks(ypos)
ax.set_yticklabels([lab for lab, _, _ in LADDER], fontsize=8.5)
ax.set_title("B   Significant region counts", loc="left", fontsize=11, fontweight="bold")
ax.legend(fontsize=8.5, frameon=False, loc="lower right")
ax.grid(axis="x", lw=0.5, alpha=0.35); ax.set_axisbelow(True)
fig.tight_layout()
rows_b = [dict(panel="B_significant_region_counts", tissue=t, scope="region", metric=lab,
               count=int(res.loc[t, k])) for lab, k, _ in LADDER for t in ("liver", "blood")]
save(fig, "fig06_dmr_overview", pd.DataFrame(rows_a + rows_b))

# ---- chromosome distribution ----
sz = pd.read_csv(os.path.join(REF_DIR, "mm39.chrom.sizes"), sep="\t", header=None,
                 names=["chrom", "len"])
sz = sz[sz.chrom.str.startswith("chr")].copy()
sz["chrom"] = sz.chrom.replace({"chrMT": "chrM"})
sz = sz.set_index("chrom")
order = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"]
order = [c for c in order if c in sz.index]

ch = {t: pd.read_csv(os.path.join(TAB, f"distribution_chrom_{t}.csv")).set_index("chrom")
      for t in ("liver", "blood")}

fig, axes = plt.subplots(2, 1, figsize=(12.6, 7.8), gridspec_kw={"height_ratios": [1.35, 1]})
rows_c = []
x = np.arange(len(order))

# --- A: raw candidate-region counts per chromosome
ax = axes[0]
for i, t in enumerate(("liver", "blood")):
    v = np.array([ch[t].reindex(order)["all_regions"].fillna(0).loc[c] for c in order])
    ax.bar(x + (0.5 - i) * 0.38, v, width=0.38, color=COL_TISSUE[t], label=t)
    for c, vv in zip(order, v):
        rows_c.append(dict(panel="A_count_per_chromosome", tissue=t, chrom=c,
                           metric="candidate_regions", value=float(vv)))
ax.set_ylabel("candidate regions")
ax.set_title("A   Candidate regions per chromosome (absolute)", loc="left", fontsize=11,
             fontweight="bold")
ax.legend(fontsize=8.5, frameon=False, loc="upper right")
ax.set_ylim(0, 12000)
ax.grid(axis="y", lw=0.5, alpha=0.35); ax.set_axisbelow(True)

# --- B: length-normalised density (chrM excluded: 16 kb, density not comparable)
order_b = [c for c in order if c != "chrM"]
ax = axes[1]
xb = np.arange(len(order_b))
for i, t in enumerate(("liver", "blood")):
    v = np.array([ch[t].reindex(order_b)["all_regions"].fillna(0).loc[c] / (sz.loc[c, "len"] / 1e6)
                  for c in order_b])
    ax.bar(xb + (0.5 - i) * 0.38, v, width=0.38, color=COL_TISSUE[t], label=t)
    for c, vv in zip(order_b, v):
        rows_c.append(dict(panel="B_density_per_Mb", tissue=t, chrom=c,
                           metric="candidate_regions_per_Mb", value=float(vv)))
ax.set_ylabel("candidate regions per Mb")
ax.set_title("B   Length-normalised density (chrM excluded: 16 kb, not comparable)",
             loc="left", fontsize=11, fontweight="bold")
ax.grid(axis="y", lw=0.5, alpha=0.35); ax.set_axisbelow(True)

# mark the single significant region (blood, chr13) with a band, no leader line
if "chr13" in order:
    xi = order.index("chr13")
    axes[0].axvspan(xi - 0.5, xi + 0.5, color="#D55E00", alpha=0.11, zorder=0)
    xbi = order_b.index("chr13")
    axes[1].axvspan(xbi - 0.5, xbi + 0.5, color="#D55E00", alpha=0.11, zorder=0)
    axes[0].text(xi, 11500, "chr13\n1 region q<0.05", ha="center", va="top",
                 fontsize=7.2, color="#D55E00", linespacing=1.25)

for ax, oo in ((axes[0], order), (axes[1], order_b)):
    ax.set_xticks(np.arange(len(oo)))
    ax.set_xticklabels([c.replace("chr", "") for c in oo], fontsize=8)
axes[0].set_xlim(-0.7, len(order) - 0.3)
axes[1].set_xlim(-0.7, len(order_b) - 0.3)
axes[1].set_xlabel("chromosome (mm39)")
fig.tight_layout()
save(fig, "fig05_chromosome_distribution", pd.DataFrame(rows_c))

print("done")


# ===========================================================================
# case-control methylation at selected loci
# ===========================================================================
FIGD = FIG                                     # figure directory
FIGNAME = {"Fig5_combined": "fig11_fig5_combined",
           "Fig5A_global_methylation": "fig08_fig5a_global_methylation",
           "Fig5B_blood_loci": "fig09_fig5b_blood_loci",
           "Fig5C_liver_loci": "fig10_fig5c_liver_loci"}

SANS = _style.apply(7.0)
plt.rcParams.update({
    "axes.linewidth": 0.7, "axes.labelsize": 7.2,
    "xtick.labelsize": 6.4, "ytick.labelsize": 6.4,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 600,
})
CB = {"Control": "#0072B2", "FLI 2h": "#D55E00"}
LBL = {"Control": "Control", "FLI 2h": "FLI 2h"}

COLS = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
        "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat",
        "cdsEndStat", "exonFrames"]
RG = pd.read_csv(os.path.join(REF, "refGene.txt.gz"), sep="\t", header=None, names=COLS,
                 low_memory=False)
CGI = pd.read_csv(os.path.join(REF, "cpgIslandExt.txt.gz"), sep="\t", header=None,
                  names=["bin", "chrom", "s", "e", "name", "length", "cpgNum", "gcNum",
                         "perCpg", "perGc", "obsExp"])

SAMP = {
    "blood": {"Control": ["Control-Blood-1", "Control-Blood-2", "Control-Blood-3"],
              "FLI 2h":  ["FAR-Blood-1", "FAR-Blood-2", "FAR-Blood-3"]},
    "liver": {"Control": ["Control-liver-1", "Control-liver-2", "Control-liver-3"],
              "FLI 2h":  ["FAR-liver-1", "FAR-liver-2", "FAR-liver-3"]},
}
LOCI = {"blood": ["Tet1", "Dnmt1", "Runx1", "Dnmt3a", "Msh2"],
        "liver": ["Dnmt1", "Hnf4a", "Fgf21"]}
COMP = {t: pd.read_csv(os.path.join(TAB, f"comparable_{t}.tsv.gz"), sep="\t") for t in SAMP}


def star(p):
    if not np.isfinite(p):
        return "n.d."
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"


def canonical_tx(gene, chrom):
    s = RG[(RG["name2"] == gene) & (RG["chrom"] == chrom)].copy()
    s["L"] = s["txEnd"] - s["txStart"]
    return s.sort_values("L").iloc[-1]


def window_for_tissue(gene, tissue):
    s = RG[RG["name2"] == gene]
    chrom = s["chrom"].iloc[0]
    tx = canonical_tx(gene, chrom)
    body_s, body_e = int(tx["txStart"]), int(tx["txEnd"])
    comp = COMP[tissue]
    for flank in (1000, 3000, 5000):
        ws, we = body_s - flank, body_e + flank
        m = (comp["chrom"] == chrom) & (comp["start"] >= ws) & (comp["end"] <= we)
        if int(m.sum()) >= 10 or flank == 5000:
            return chrom, ws, we, int(m.sum()), flank, tx
    raise RuntimeError


def locus_data(gene, tissue):
    chrom, ws, we, n_comp, flank, tx = window_for_tissue(gene, tissue)
    comp = COMP[tissue]
    m = (comp["chrom"] == chrom) & (comp["start"] >= ws) & (comp["end"] <= we)
    sel = comp[m].sort_values("start")
    pos = sel["start"].to_numpy()
    pooled, per = {}, {}
    for g, sids in SAMP[tissue].items():
        M = sel[[x + "_M" for x in sids]].to_numpy(float)
        C = sel[[x + "_cov" for x in sids]].to_numpy(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            pooled[g] = (np.where(C.sum(1) > 0, M.sum(1) / np.where(C.sum(1) > 0, C.sum(1), 1), np.nan),
                         C.sum(1))
            vals = np.where(C > 0, M / np.where(C > 0, C, 1), np.nan)
        per[g] = np.nanmean(vals, axis=0)
    c, f = per["Control"], per["FLI 2h"]
    ok = np.isfinite(c).all() and np.isfinite(f).all() and len(c) == 3 and len(f) == 3
    p = stats.ttest_ind(f, c, equal_var=False).pvalue if ok else np.nan
    rec = dict(tissue=tissue, gene=gene, chrom=chrom, window_start=ws, window_end=we,
               window_kb=round((we - ws) / 1000, 2), flank_bp=flank,
               n_comparable_cpg=n_comp,
               control_mean=float(np.nanmean(c)), fli2h_mean=float(np.nanmean(f)),
               delta=float(np.nanmean(f) - np.nanmean(c)),
               control_sd=float(np.nanstd(c, ddof=1)), fli2h_sd=float(np.nanstd(f, ddof=1)),
               welch_p=float(p) if np.isfinite(p) else np.nan, mark=star(p),
               control_samples=";".join(f"{v:.4f}" for v in c),
               fli2h_samples=";".join(f"{v:.4f}" for v in f))
    return sel, pos, pooled, rec, tx


# ---------------- primitives ----------------
def pos_fmt(ws, we):
    """span-adaptive coordinate format: kb below 1 Mb, Mb above (avoids duplicate labels)"""
    span = we - ws
    if span < 1_000_000:
        return FuncFormatter(lambda x, _: f"{x/1e3:,.0f}"), "kb"
    return FuncFormatter(lambda x, _: f"{x/1e6:.2f}"), "Mb"


def draw_track(ax, pos, beta, cov, color, ws, we, label=None, show_ticks=True):
    keep = np.isfinite(beta) & (cov > 0)
    ax.vlines(pos[keep], 0, beta[keep], color=color, lw=1.25, alpha=0.92, zorder=3)
    ax.set_xlim(ws, we); ax.set_ylim(0, 1)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"] if show_ticks else ["", ""])
    ax.tick_params(axis="both", length=2, pad=1.2)
    ax.tick_params(axis="x", labelbottom=False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.5); ax.spines[sp].set_color("#999999")
    if label:
        ax.set_ylabel(label, fontsize=6.2, rotation=0, ha="right", va="center", labelpad=8,
                      color=color)


def draw_genemodel(ax, tx, ws, we):
    a_s = np.fromstring(tx["exonStarts"].rstrip(","), sep=",", dtype=np.int64)
    a_e = np.fromstring(tx["exonEnds"].rstrip(","), sep=",", dtype=np.int64)
    ax.hlines(0.55, max(ws, int(tx["txStart"])), min(we, int(tx["txEnd"])),
              color="#333333", lw=0.9, zorder=2)
    minw = (we - ws) * 0.0022
    for a, b in zip(a_s, a_e):
        if b < ws or a > we:
            continue
        ax.add_patch(Rectangle((max(a, ws), 0.18), max(b - a, minw), 0.74,
                               color="#333333", lw=0, zorder=3))
    strand = tx["strand"]; d = 1 if strand == "+" else -1
    tip = min(we, int(tx["txEnd"])) if d == 1 else max(ws, int(tx["txStart"]))
    ax.add_patch(FancyArrow(tip - d * (we - ws) * 0.010, 0.55, d * (we - ws) * 0.010, 0,
                            width=(we - ws) * 0.0009, head_width=(we - ws) * 0.0055,
                            head_length=(we - ws) * 0.007, color="#333333", lw=0,
                            length_includes_head=True, zorder=4))
    ax.set_xlim(ws, we); ax.set_ylim(0, 1); ax.set_yticks([])
    ax.tick_params(axis="x", length=2, labelsize=6.3, pad=1.2)
    f, unit = pos_fmt(ws, we)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    ax.xaxis.set_major_formatter(f)
    ax.set_xlabel(unit, fontsize=5.6, labelpad=1.0)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.5); ax.spines[sp].set_color("#999999")


def draw_cgi(ax, chrom, ws, we):
    sub = CGI[(CGI["chrom"] == chrom) & (CGI["e"] > ws) & (CGI["s"] < we)]
    for _, r in sub.iterrows():
        ax.add_patch(Rectangle((max(r["s"], ws), 0.05), min(r["e"], we) - max(r["s"], ws), 0.9,
                               color="#009E73", lw=0, alpha=0.8, zorder=2))
    ax.set_xlim(ws, we); ax.set_ylim(0, 1); ax.set_yticks([]); ax.set_xticks([])
    for sp in ("left",):
        ax.spines[sp].set_linewidth(0.5); ax.spines[sp].set_color("#999999")
    ax.spines["bottom"].set_visible(False)


def draw_minibar(ax, rec):
    vals = [rec["control_mean"], rec["fli2h_mean"]]
    sds = [rec["control_sd"], rec["fli2h_sd"]]
    xs = [0, 1]
    ax.bar(xs, vals, width=0.60, color=[CB["Control"], CB["FLI 2h"]], edgecolor="none")
    ax.errorbar(xs, vals, yerr=sds, fmt="none", ecolor="#333333", elinewidth=0.7,
                capsize=2.0, capthick=0.7)
    ax.set_ylim(0, 1.18); ax.set_xlim(-0.65, 1.65)
    ax.set_xticks(xs); ax.set_xticklabels(["C", "F"], fontsize=6.2)
    ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(["0", "0.5", "1"], fontsize=6.2)
    ax.tick_params(length=2, pad=1.0)
    top = max(np.array(vals) + np.array(sds))
    y = min(top + 0.09, 0.96)
    ax.plot([0, 0, 1, 1], [y, y + 0.04, y + 0.04, y], color="#333333", lw=0.6)
    ax.text(0.5, y + 0.055, rec["mark"], ha="center", va="bottom", fontsize=7.2)
    ax.set_title("p=%s" % (f"{rec['welch_p']:.2g}" if np.isfinite(rec["welch_p"]) else "n.d."),
                 fontsize=5.6, pad=1.5, color="#444444")
    ax.spines["right"].set_visible(False)


def txt(fig, W, H, x_in, y_in, s, **kw):
    """text placed with inch coordinates measured from the top-left corner"""
    fig.text(x_in / W, 1 - y_in / H, s, **kw)


def ax_in(fig, W, H, x, y, w, h):
    """axes rect from inch coordinates measured from the top-left corner"""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


# ---------------- panels ----------------
def panel_A(fig, W, H, x, y):
    gq = pd.read_csv(os.path.join(TAB, "sample_QC.csv"))
    gm = pd.read_csv(os.path.join(TAB, "global_methylation.csv"))
    specs = [("A1", "all CpGs", gq, "global_beta_unweighted"),
             ("A2", "CpGs covered ≥5× in ≥2/3 per group", gm, "mean_beta_unweighted")]
    rows, prow = [], []
    for k, (tag, label, df, col) in enumerate(specs):
        ax = ax_in(fig, W, H, x + k * 2.32, y, 2.06, 1.50)
        xs, hs, es, xt, xl = [], [], [], [], []
        for i, tis in enumerate(("blood", "liver")):
            for j, g in enumerate(("Control", "FLI 2h")):
                sel = df[(df["tissue"] == tis) &
                         (df["group"] == ("Control" if g == "Control" else "FAR"))]
                v = sel[col].to_numpy(float)
                xb = i * 2.2 + j * 0.85
                xs.append(xb); hs.append(v.mean() * 100); es.append(v.std(ddof=1) * 100)
                rows += [dict(panel="A", definition=tag, metric=label, tissue=tis,
                              group=g, sample_id=sid, beta_fraction=val)
                         for sid, val in zip(sel["sample_id"], v)]
            c = df[(df["tissue"] == tis) & (df["group"] == "Control")][col].to_numpy(float)
            f = df[(df["tissue"] == tis) & (df["group"] == "FAR")][col].to_numpy(float)
            p = stats.ttest_ind(f, c, equal_var=False).pvalue
            prow.append(dict(panel="A", definition=tag, tissue=tis,
                             control_mean_pct=c.mean() * 100, fli2h_mean_pct=f.mean() * 100,
                             delta_pct=(f.mean() - c.mean()) * 100,
                             welch_p=p, mark=star(p), n_per_group=len(c)))
            xt.append(i * 2.2 + 0.425); xl.append(tis.capitalize())
            yl = max(np.array(hs) + np.array(es)) + 6
            ax.plot([xt[-1] - 0.36, xt[-1] - 0.36, xt[-1] + 0.36, xt[-1] + 0.36],
                    [yl, yl + 3, yl + 3, yl], color="#333333", lw=0.6)
            ax.text(xt[-1], yl + 4, f"{star(p)}  (p={p:.2f})", ha="center", va="bottom",
                    fontsize=6.0)
        ax.bar(xs, hs, width=0.72,
               color=[CB["Control"] if i % 2 == 0 else CB["FLI 2h"] for i in range(len(xs))],
               edgecolor="none")
        ax.errorbar(xs, hs, yerr=es, fmt="none", ecolor="#333333", elinewidth=0.7,
                    capsize=2.2, capthick=0.7)
        ax.set_xticks(xt); ax.set_xticklabels(xl)
        ax.set_ylabel("CpG methylation (%)"); ax.set_ylim(0, 100)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_title(f"{tag}   {label}", loc="left", fontsize=6.8, fontweight="bold", pad=3)
        if k == 1:
            ax.set_ylabel("")
            txt(fig, W, H, 5.28, y - 0.02, "", fontsize=6.0)
            fig.legend(handles=[Rectangle((0, 0), 1, 1, color=CB["Control"], label=LBL["Control"]),
                                Rectangle((0, 0), 1, 1, color=CB["FLI 2h"], label=LBL["FLI 2h"])],
                       fontsize=6.2, frameon=False, ncol=1, loc="upper left",
                       bbox_to_anchor=(5.22 / W, 1 - (y - 0.04) / H), handlelength=1.1,
                       borderpad=0.1, labelspacing=0.3)
        ax.grid(axis="y", lw=0.4, alpha=0.3); ax.set_axisbelow(True)
    return pd.DataFrame(rows), pd.DataFrame(prow)


def panel_loci(fig, W, H, tissue, genes, x, y, label):
    w_track, h_t1, h_t2, h_model, h_cgi = 4.52, 0.30, 0.30, 0.15, 0.07
    recs, rows = [], []
    for k, gene in enumerate(genes):
        sel, pos, pooled, rec, tx = locus_data(gene, tissue)
        ws, we, chrom = rec["window_start"], rec["window_end"], rec["chrom"]
        ax1 = ax_in(fig, W, H, x, y, w_track, h_t1)
        draw_track(ax1, pos, pooled["Control"][0], pooled["Control"][1], CB["Control"], ws, we,
                   label=LBL["Control"], show_ticks=(k == 0))
        # chromosome, window length and CpG count; no exact coordinates, no qualifier
        ax1.set_title(f"{gene}   {chrom}   ({rec['window_kb']:.1f} kb, "
                      f"{rec['n_comparable_cpg']} CpGs)",
                      loc="left", fontsize=6.6, fontweight="bold", pad=2.2)
        ax2 = ax_in(fig, W, H, x, y + h_t1 + 0.035, w_track, h_t2)
        draw_track(ax2, pos, pooled["FLI 2h"][0], pooled["FLI 2h"][1], CB["FLI 2h"], ws, we,
                   label=LBL["FLI 2h"], show_ticks=(k == 0))
        # CpG strip ABOVE the gene model: drawn below it, the strip overpaints the band
        # where the coordinate tick labels live and clips their glyphs
        axc = ax_in(fig, W, H, x, y + h_t1 + h_t2 + 0.035 + 0.045, w_track, h_cgi)
        draw_cgi(axc, chrom, ws, we)
        axm = ax_in(fig, W, H, x, y + h_t1 + h_t2 + 0.035 + 0.045 + h_cgi, w_track, h_model)
        draw_genemodel(axm, tx, ws, we)
        axb = ax_in(fig, W, H, x + w_track + 0.34, y, 0.60, h_t1 + h_t2 + 0.035)
        draw_minibar(axb, rec)
        if k == 0:
            axb.set_ylabel("methylation", fontsize=6.0, labelpad=2)
            axc.text(-0.012, 0.5, "CpG island", transform=axc.transAxes, fontsize=5.2,
                     color="#009E73", ha="right", va="center")
        recs.append(rec)
        for g in ("Control", "FLI 2h"):
            for sid, v in zip(SAMP[tissue][g],
                              [float(t) for t in (rec["control_samples"] if g == "Control"
                                                  else rec["fli2h_samples"]).split(";")]):
                rows.append(dict(panel=label, tissue=tissue, gene=gene, chrom=chrom,
                                 window_start=ws, window_end=we, n_comparable_cpg=rec["n_comparable_cpg"],
                                 group=g, sample_id=sid, window_mean_beta=v))
        y += h_t1 + h_t2 + 0.035 + 0.045 + h_model + h_cgi + 0.30
    return recs, rows, y


def build(panels, name):
    H = 0.30
    if "A" in panels: H += 1.50 + 0.42
    if "B" in panels: H += 0.30 + 5 * 1.17
    if "C" in panels: H += 0.30 + 3 * 1.17
    H += 0.82
    W = 6.95
    fig = plt.figure(figsize=(W, H))
    y = 0.30
    x0 = 0.72
    recA, rowA, recL, rowL = [], [], [], []
    if "A" in panels:
        txt(fig, W, H, 0.30, y + 0.16, "A", fontsize=11.5, fontweight="bold")
        ra, pa = panel_A(fig, W, H, x0, y)
        recA.append(pa); rowA.append(ra)
        y += 1.50 + 0.42
    if "B" in panels:
        txt(fig, W, H, 0.10, y + 0.05, "B", fontsize=11.5, fontweight="bold")
        txt(fig, W, H, 0.36, y + 0.055, "Blood", fontsize=8.8, fontweight="bold", color="#333333")
        fig.add_artist(plt.Line2D([0.72 / W, 5.24 / W], [1 - (y + 0.115) / H] * 2,
                                  color="#bbbbbb", lw=0.7, transform=fig.transFigure))
        r, w, y = panel_loci(fig, W, H, "blood", LOCI["blood"], x0, y + 0.28, "B")
        recL += r; rowL += w
        y += 0.10
    if "C" in panels:
        txt(fig, W, H, 0.10, y + 0.05, "C", fontsize=11.5, fontweight="bold")
        txt(fig, W, H, 0.36, y + 0.055, "Liver", fontsize=8.8, fontweight="bold", color="#333333")
        fig.add_artist(plt.Line2D([0.72 / W, 5.24 / W], [1 - (y + 0.115) / H] * 2,
                                  color="#bbbbbb", lw=0.7, transform=fig.transFigure))
        r, w, y = panel_loci(fig, W, H, "liver", LOCI["liver"], x0, y + 0.28, "C")
        recL += r; rowL += w
    foot = ("Bars: group mean ± SD over the three animals (n = 3 vs 3), Welch t-test (nominal); "
            "exact permutation p cannot reach 0.05 at 3 vs 3.")
    if "B" in panels or "C" in panels:
        foot = ("Tracks: merged group signal - per-CpG methylation of the three animals pooled by "
                "reads (0-1); green strip = CpG island (absent = none annotated in that window). " + foot)
    txt(fig, W, H, 0.42, H - 0.34, foot, fontsize=5.6, color="#444444", wrap=True)
    pdf = os.path.join(FIGD, FIGNAME[name] + ".pdf")
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", pdf, f"({W:.2f}x{H:.2f} in)")
    return recA, rowA, recL, rowL


if __name__ == "__main__":
    A, RA, L, RL = [], [], [], []
    for panels, name in ((("A", "B", "C"), "Fig5_combined"),
                         (("A",), "Fig5A_global_methylation"),
                         (("B",), "Fig5B_blood_loci"),
                         (("C",), "Fig5C_liver_loci")):
        a, ra, l, rl = build(panels, name)
        A += a; RA += ra; L += l; RL += rl
    st = pd.DataFrame(L).drop_duplicates(subset=["tissue", "gene"])
    st.to_csv(os.path.join(TAB, "fig5_locus_statistics.csv"), index=False, float_format="%.6g")
    sa = pd.concat(A, ignore_index=True).drop_duplicates(subset=["definition", "tissue"])
    sa.to_csv(os.path.join(TAB, "fig5a_global_statistics.csv"), index=False, float_format="%.6g")
    pd.DataFrame(RL).drop_duplicates(["tissue", "gene", "group", "sample_id"]) \
        .to_csv(os.path.join(TAB, "fig5_locus_samples.csv"), index=False, float_format="%.6g")
    pd.concat(RA, ignore_index=True).drop_duplicates(["definition", "tissue", "group", "sample_id"]) \
        .to_csv(os.path.join(TAB, "fig5a_global_samples.csv"), index=False, float_format="%.6g")
    print("\n--- panel A ---")
    print(sa[["definition", "tissue", "control_mean_pct", "fli2h_mean_pct", "delta_pct",
              "welch_p", "mark"]].to_string(index=False))
    print("\n--- loci ---")
    print(st[["tissue", "gene", "chrom", "window_kb", "n_comparable_cpg", "control_mean",
              "fli2h_mean", "delta", "welch_p", "mark"]].to_string(index=False))


# ===========================================================================
# locus panels: epigenetic regulators and functional genes
# ===========================================================================
_style.apply(7.0)

TISSUE = "blood"

PANELS = [
    ("fig12_supp_fig4a_epigenetic_loci", "A", ["Dnmt3a", "Dnmt3b", "Tet1", "Tet2"]),
    ("fig13_supp_fig4b_functional_loci", "B", ["Fgf21", "Hnf4a", "Msh2", "Runx1", "Trp53"]),
]

FOOT = ("Tracks: merged group signal - per-CpG methylation of the three animals pooled by reads "
        "(0-1); green strip = CpG island (absent = none annotated in that window). "
        "Bars: group mean \u00b1 SD over the three animals (n = 3 vs 3), Welch t-test (nominal); "
        "exact permutation p cannot reach 0.05 at 3 vs 3.")

all_rec, all_rows = [], []
for fname, label, genes in PANELS:
    W = 6.95
    H = 0.30 + 0.30 + len(genes) * 1.17 + 0.82
    fig = plt.figure(figsize=(W, H))
    y = 0.30
    x0 = 0.72
    txt(fig, W, H, 0.10, y + 0.05, label, fontsize=11.5, fontweight="bold")
    txt(fig, W, H, 0.36, y + 0.055, "Blood", fontsize=8.8, fontweight="bold", color="#333333")
    fig.add_artist(plt.Line2D([x0 / W, 5.24 / W], [1 - (y + 0.115) / H] * 2,
                              color="#bbbbbb", lw=0.7, transform=fig.transFigure))

    recs, rows, y_end = panel_loci(fig, W, H, TISSUE, genes, x0, y + 0.28, label)
    txt(fig, W, H, 0.42, H - 0.34, FOOT, fontsize=5.6, color="#444444", wrap=True)

    pdf = os.path.join(FIG, fname + ".pdf")
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", pdf, f"({W:.2f}x{H:.2f} in = {W*25.4:.0f}x{H*25.4:.0f} mm)")

    d = pd.DataFrame(rows)
    d.insert(0, "figure", fname)
    d.drop_duplicates(["figure", "gene", "group", "sample_id"]) \
        .to_csv(os.path.join(FIG, fname + ".csv"), index=False, float_format="%.6g")

    for r in recs:
        r = dict(r); r["figure"] = fname; all_rec.append(r)

st = pd.DataFrame(all_rec)
st.to_csv(os.path.join(TAB, "supp_fig4_locus_statistics.csv"), index=False, float_format="%.6g")
print()
print(st[["figure", "gene", "chrom", "window_kb", "n_comparable_cpg", "control_mean",
          "fli2h_mean", "delta", "welch_p", "mark"]].to_string(index=False))


# ===========================================================================
# imprinting-control-region panel
# ===========================================================================
_style.apply(7.0)

GENES_ALL = ["H19", "Kcnq1ot1", "Peg3", "Ndn", "Impact", "Gnas"]
GENES = [g for g in GENES_ALL if g != "Impact"]   # Impact is computed but not displayed
PANELS = [("fig14_imprinted_loci_blood", "A", "blood"),
          ("fig15_imprinted_loci_liver", "B", "liver")]

UP, DN = 6000, 2000
MIN_CPG = 10
WIDEN = [(UP, DN), (10000, 4000), (20000, 8000)]


def tss_window(gene, tissue):
    """TSS-anchored window [-UP, +DN] of the longest transcript, widened for coverage."""
    s = RG[RG.name2 == gene]
    chrom = s["chrom"].iloc[0]
    tx = canonical_tx(gene, chrom)
    tss = int(tx["txStart"]) if tx["strand"] == "+" else int(tx["txEnd"])
    comp = COMP[tissue]
    for up, dn in WIDEN:
        ws, we = (tss - up, tss + dn) if tx["strand"] == "+" else (tss - dn, tss + up)
        m = (comp["chrom"] == chrom) & (comp["start"] >= ws) & (comp["end"] <= we)
        if int(m.sum()) >= MIN_CPG or (up, dn) == WIDEN[-1]:
            return chrom, ws, we, int(m.sum()), up, tx
    raise RuntimeError


window_for_tissue = tss_window          # panel_loci / locus_data pick this up

FOOT = ("Tracks: merged group signal - per-CpG methylation of the three animals pooled by reads "
        "(0-1); green strip = annotated CpG island. Windows are TSS-anchored (-6 kb / +2 kb) "
        "proxies for the promoter-proximal germline DMRs; element-level ICR boundaries were not "
        "available offline. Bars: group mean \u00b1 SD (n = 3 vs 3), Welch t-test. P values are "
        "nominal and uncorrected for the 12 windows examined (exact permutation floor P = 0.10 at "
        "3 vs 3); the overlapping candidate region(s) reach region-level q = 0.97-1.00 in every "
        "window (table), so no ICR window is significant under the prespecified BH framework.")

all_rec, all_rows = [], []
for fname, label, tissue in PANELS:
    W = 6.95
    H = 0.30 + 0.30 + len(GENES) * 1.17 + 0.82
    fig = plt.figure(figsize=(W, H))
    y, x0 = 0.30, 0.72
    txt(fig, W, H, 0.10, y + 0.05, label, fontsize=11.5, fontweight="bold")
    txt(fig, W, H, 0.36, y + 0.055, tissue.capitalize(), fontsize=8.8, fontweight="bold",
          color="#333333")
    fig.add_artist(plt.Line2D([x0 / W, 5.24 / W], [1 - (y + 0.115) / H] * 2,
                              color="#bbbbbb", lw=0.7, transform=fig.transFigure))

    recs, rows, y_end = panel_loci(fig, W, H, tissue, GENES, x0, y + 0.28, label)
    txt(fig, W, H, 0.42, H - 0.40, FOOT, fontsize=5.5, color="#444444", wrap=True)

    pdf = os.path.join(FIG, fname + ".pdf")
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", pdf, f"({W:.2f}x{H:.2f} in = {W*25.4:.0f}x{H*25.4:.0f} mm)")

    d = pd.DataFrame(rows)
    d.insert(0, "figure", fname)
    d.drop_duplicates(["figure", "gene", "group", "sample_id"]) \
        .to_csv(os.path.join(FIG, fname + ".csv"), index=False, float_format="%.6g")

    for r in recs:
        r = dict(r)
        r["figure"] = fname
        r["shown_in_figure"] = True
        all_rec.append(r)
    # loci computed for the record but deliberately not drawn
    for g in [x for x in GENES_ALL if x not in GENES]:
        _, _, _, rec, _ = locus_data(g, tissue)
        rec = dict(rec)
        rec["figure"] = fname
        rec["shown_in_figure"] = False
        all_rec.append(rec)

st = pd.DataFrame(all_rec)

# ---- region-level evidence for the same windows (kept separate from the Welch test) ----
rows_q = []
for tissue in sorted(st["tissue"].unique()):
    reg = pd.read_csv(os.path.join(TAB, f"regions_all_{tissue}.tsv.gz"), sep="\t")
    tst = pd.read_csv(os.path.join(TAB, f"region_tests_{tissue}.tsv.gz"), sep="\t")
    qcol = [c for c in tst.columns if c.startswith("q_primary")][0]
    for _, r in st[st.tissue == tissue].iterrows():
        ov = reg[(reg.chrom == r.chrom) & (reg.end > r.window_start) & (reg.start < r.window_end)]
        q = ov.merge(tst[["region_id", qcol]], on="region_id", how="left")
        rows_q.append(dict(figure=r.figure, tissue=tissue, gene=r.gene,
                           n_candidate_regions=len(ov),
                           min_region_q=(float(q[qcol].min()) if len(q) else float("nan")),
                           region_logFC_max=(float(ov["diff"].abs().max()) if len(ov) else float("nan"))))
st = st.merge(pd.DataFrame(rows_q), on=["figure", "tissue", "gene"], how="left")
st.to_csv(os.path.join(TAB, "icr_locus_statistics.csv"), index=False, float_format="%.6g")
print()
print(st[["figure", "gene", "shown_in_figure", "chrom", "window_kb", "n_comparable_cpg",
          "control_mean", "fli2h_mean", "delta", "welch_p", "n_candidate_regions",
          "min_region_q"]].to_string(index=False))


print("\nall methylation figures written to", FIG)
