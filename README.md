# Methylation

DNA methylation of offspring blood and liver, FLI 2 h versus unexposed control,
mouse mm39. Blood and liver are analysed separately, three samples per group.
The effect size is FLI 2 h minus control, so a positive value means higher
methylation in the imaged group.

## Scripts, in run order

| # | Script | Purpose |
|---|---|---|
| 1 | `cpg_counts_qc.py` | per-sample CpG tables → coverage QC, comparable CpG sets, global methylation |
| 2 | `verify_cpg_columns.py` | independent check of the count columns and coordinates |
| 3 | `build_regions.py` | candidate regions from coverage only; per-sample region methylation and effect size |
| 4 | `region_test_limma.R` | region-level test, primary and sensitivity models |
| 5 | `dss_dml.R` | DSS site-level and region-level branch |
| 6 | `annotate_regions.py` | genomic annotation and result layers |
| 7 | `supplementary_checks.py` | coverage sensitivity, DSS versus region-test agreement, locus re-check |
| 8 | `figures.py` | every figure, fig01-fig15: sample QC and effect sizes, region and chromosome distribution, Fig. 5A-5C, Supplementary Fig. 4A-4B, ICR panel |

`plot_style.py` is a helper imported by `figures.py`.

## Data locations

Each script starts with a block to edit:

```python
IN_CPG_DIR = "data/cpg_tables"   # per-sample CpG tables
REF_DIR    = "data/mm39"         # mm39 annotation
OUT_DIR    = "results"           # tables/ and figures/
```

| Input | Format |
|---|---|
| Per-sample CpG tables | one file per sample, no header, tab-separated: `chrom, start0, end, beta, M, U` |
| mm39 annotation | `mm39.chrom.sizes`, `refGene.txt.gz`, `cpgIslandExt.txt.gz` (UCSC) |

To start from the deposited processed tables instead of step 1, put
`comparable_liver.tsv.gz`, `comparable_blood.tsv.gz` and the step-1 QC tables in
`results/tables/` and begin at step 3.

## Statistical units

The DSS branch, the region-level test and the coverage-defined candidate regions are
three different units. Their counts are reported separately and are never added
together. Merged group tracks in the figures are for display; every statistic uses the
individual samples.

## Reproduction start

Existing per-sample CpG counts. Alignment and methylation calling are not part of this
branch.

## Software

Python 3.10 with numpy 2.2.6, pandas 2.3.3, scipy 1.15.3, matplotlib 3.8.4.
R 4.5.3 with DSS 2.58.0, bsseq 1.46.0, GenomicRanges 1.62.1, limma 3.66.0,
data.table 1.18.6.1. Reference: UCSC mm39.
