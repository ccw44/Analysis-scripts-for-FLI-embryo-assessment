# Mouse RNA-seq preprocessing

Upstream scripts will be adapted from the [RNAseq branch of the previous BMAL1 study](https://github.com/ccw44/Analysis-scripts-for-BMAL1-mediated-neurodifferentiation/tree/RNAseq).

## Workflow

Read QC → adapter/quality trimming → genome alignment → BAM sorting/indexing → gene-level counts.

Use mouse GRCm39/mm39 genome and matching annotation resources instead of the human references in the original scripts. Source scripts cover FastQC/MultiQC, Trim Galore, HISAT2 or STAR, SAMtools and featureCounts. The study manuscript currently specifies HTSeq. The actual alignment and counting tools must be confirmed from the original analysis records before the adapted code is finalized; changing the README alone does not establish which workflow produced the published results.

## Scripts pending preparation

Adapt the necessary QC, trimming, alignment, BAM and counting scripts. Retain TPM/FPKM conversion only if required by the final analysis. Remove private paths and unrelated study examples. Preserve actual analysis settings in the scripts and document the software versions.

The mouse genome assembly has been confirmed by the author; the annotation release, selected alignment workflow and counting tool remain to be verified. Adaptation and data-level reproduction have not yet been performed. Downstream analysis scripts from the previous study are not automatically adopted.
