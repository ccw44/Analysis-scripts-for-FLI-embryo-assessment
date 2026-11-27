# Proteomics

Protein abundance in high-score (HS) versus low-score (LS) embryos, independent
two-group design, three embryos per group. A positive log2 fold change means higher
abundance in HS.

## Scripts, in run order

| # | Script | Purpose |
|---|---|---|
| 1 | `prepare_inputs.py` | stage the vendor package, write the manifest and sample sheet |
| 2 | `analysis.R` | QC, filtering, log2 transform, limma moderated test, BH |
| 3 | `enrichment.R` | gene collapse, GO / KEGG over-representation and GSEA |
| 4 | `figures.R` | every figure: Fig. 2B, Fig. 2C, the GSEA and per-protein supplements, the supplementary ORA panels, and the two exploratory panels |

```bash
python prepare_inputs.py results/proteomics
Rscript analysis.R  results/proteomics
Rscript enrichment.R results/proteomics
Rscript figures.R      results/proteomics
```

## Two branches, never merged

| Branch | Steps | Own candidate set and figures |
|---|---|---|
| Matrix reanalysis | 2 → 3 | `statistics/main_all.tsv`, `enrichment/*`, Fig. 2B, GSEA panels |
| Provider result review | 1 | `audit/legacy_statistics_BH_recomputed.tsv`, Fig. 2C |

Both are drawn by `figures.R`, in separate sections that keep their own
candidate set and enrichment background. Two further sections are exploratory: they use
a raw P < 0.05 and |log2FC| ≥ 1 foreground, which is not a multiple-testing-supported
result and is not presented as one.

## Input

The vendor `Result_N0011` directory, placed at `<base>/inputs/Result_N0011`, or the
original archive via `--archive`. Raw mass-spectrometry files are not available; the
matrix is the starting point.

## Model

Zero means missing; log2; no imputation for the primary inference. Filter: at least two
observed values in each group. limma `lmFit` + `eBayes(robust = TRUE)`, HS minus LS, BH
FDR < 0.05, with |log2FC| ≥ 1 as a separate effect-size subset.

## Software

R 4.5.3 with limma 3.66.0, data.table 1.18.6.1, org.Mm.eg.db 3.22.0, GO.db 3.22.0,
fgsea 1.36.2, ggplot2 4.0.3, patchwork 1.3.2, ggrepel 0.9.8. Python 3.10 for step 1.
KEGG pathways are fetched from `rest.kegg.jp` at run time and are therefore not frozen;
GO terms are fixed by the pinned `GO.db`.
