# Analysis scripts for FLI embryo assessment

Bioinformatics analysis for the study **Femtosecond label-free optical sensing of embryo developmental competence using an endogenous NAD(P)H/lipofuscin ratio**.

This repository organizes RNA-seq, DNA methylation and proteomics analysis code for mouse embryos and offspring. Analysis scripts will be uploaded separately after preparation and verification.

## Analysis branches

| Branch | Analysis |
|---|---|
| [RNAseq](https://github.com/ccw44/Analysis-scripts-for-FLI-embryo-assessment/tree/RNAseq) | RNA-seq preprocessing adapted from the previous study to mouse GRCm39/mm39 |
| [Methylation](https://github.com/ccw44/Analysis-scripts-for-FLI-embryo-assessment/tree/Methylation) | Global methylation, regional differences, annotation and visualization in offspring blood and liver |
| [Proteomics](https://github.com/ccw44/Analysis-scripts-for-FLI-embryo-assessment/tree/Proteomics) | Protein-matrix processing, differential abundance, functional enrichment and visualization |

Each branch contains a brief description; analysis scripts are pending upload. This main branch contains the project overview and license.

## Data and reproducibility

RNA-seq and methylation sequencing files have been uploaded to GSA according to the study authors. Accession numbers and access links will be added after verification. Processed proteomics matrices are handled separately from this code repository. Original mass-spectrometry files have not been provided by the service company.

Methylation reanalysis starts from per-sample CpG counts; proteomics analysis starts from protein quantification matrices. No completed end-to-end reproduction is claimed for scripts that have not yet been uploaded or tested.

## Reference workflow

The RNA-seq adaptation is based on the [RNAseq branch of the previous BMAL1 study](https://github.com/ccw44/Analysis-scripts-for-BMAL1-mediated-neurodifferentiation/tree/RNAseq). The mouse reference is GRCm39/mm39, with matching genome and annotation resources. The original analysis records must resolve the featureCounts-versus-HTSeq discrepancy before the final counting script and manuscript description are aligned.

## License

Code is distributed under the MIT License. Third-party code retains its original attribution and license. The code license does not automatically apply to datasets.
