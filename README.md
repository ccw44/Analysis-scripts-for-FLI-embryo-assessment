# RNA-seq upstream processing

Raw reads to gene-level counts, mouse GRCm39. Three steps: QC, trimming, alignment
and counting.

## Scripts, in run order

| # | Script | Purpose |
|---|---|---|
| 1 | `qc.sh` | FastQC per file, then MultiQC over the set |
| 2 | `cleanfq.sh` | adapter and quality trimming with Trim Galore (paired-end) |
| 3 | `map2.sh` | HISAT2 or STAR alignment, coordinate sort, index |
| 4 | `bai.sh` | index any BAM the alignment step produced |
| 5 | `count.sh` | featureCounts gene-level and transcript-level counts |
| 6 | `count2tpm.R` | TPM / FPKM from the featureCounts output (display only) |

```bash
bash qc.sh      data/fastq            results/qc
bash cleanfq.sh data/fastq            results
bash map2.sh    results/cleanfq       results/bam  STAR
bash bai.sh     results/bam
bash count.sh   results/bam           results/count STAR
cd results/count && Rscript ../../count2tpm.R
```

`cleanfq.sh` writes `<sample>_1_val_1.fq.gz` / `<sample>_2_val_2.fq.gz`, and the later
steps read sample names back from exactly those files, so the handoff cannot drift.

## Setup

Set the reference resources before running `map2.sh` and `count.sh`:

```bash
export HISAT2_INDEX=/path/to/GRCm39/hisat2_index/Mus_musculus.GRCm39.dna.primary_assembly
export STAR_INDEX=/path/to/GRCm39/star_index
export GTF=/path/to/GRCm39/Mus_musculus.GRCm39.111.gtf
```

The paths in the scripts point at GRCm39 release 111, chosen to mirror the release
used in the reference repository (human GRCh38.111). Annotation release and
chromosome naming have to match between the index and the GTF: Ensembl names
chromosomes `1`, `2`, ... `X`, `Y`, `MT`, while UCSC-style resources use `chr1`,
`chr2`, ... A BAM built on one naming and counted with the other yields zero counts.

## Items to verify against the original analysis record

The upstream scripts were adapted from the reference repository, not rebuilt from the
manuscript text. Before these are treated as the record of what was run:

1. **Aligner actually used.** `map2.sh` supports both HISAT2 and STAR; the branch does
   not record which one produced the published counts.
2. **Counter.** `count.sh` uses featureCounts. The current manuscript text says HTSeq.
   The code was not changed to match the text; one of the two has to be corrected.
3. **Annotation release.** Release 111 is an assumption from the reference repository.
4. **Chromosome naming** as described above.
5. **`featureCounts -p`.** From subread 2.0.2, `-p` alone no longer counts fragments;
   `--countReadPairs` is required for that meaning. The flag was left as it stands and
   the subread version should be checked before comparing counts.
6. **Strandedness.** No `-s` value is set in `count.sh`; confirm the library protocol.

## Status

Syntax-checked only. These scripts have not been run on this study's data in this
workspace: no FASTQ, BAM or count file for the RNA-seq experiment is present here, and
no resource index is installed. They are not a substitute for the original run record.

## Software

FastQC, MultiQC, Trim Galore, HISAT2 or STAR, samtools, subread (featureCounts), R with
base only for `count2tpm.R`. Versions are not fixed by this branch; record the versions
of the environment that produced the counts.
