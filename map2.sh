#!/bin/bash
# Align trimmed reads to GRCm39 and write a coordinate-sorted, indexed BAM.
#
#   bash map2.sh <clean_dir> <bam_dir> <hisat2|STAR>
#
# Input   <clean_dir>/<sample>_1_val_1.fq.gz, <sample>_2_val_2.fq.gz  (from cleanfq.sh)
# Output  hisat2: <bam_dir>/<sample>.Hisat_aln.sorted.bam (+ .bai)
#         STAR:   <bam_dir>/<sample>.STAR.Aligned.sortedByCoord.out.bam (+ .bai)
#
# Sample names are read back from the trimmed files, so this step and cleanfq.sh
# cannot drift apart. Set HISAT2_INDEX / STAR_INDEX to the GRCm39 resources.
#
# Reference used for the resource paths below: GRCm39, Ensembl release 111,
# chosen to mirror the release used in the reference repository (human GRCh38.111).
# Confirm the release and the chromosome naming against the original analysis
# record before comparing with the published counts.
set -euo pipefail

clean_dir=${1:?usage: map2.sh <clean_dir> <bam_dir> <hisat2|STAR>}
bam_dir=${2:?usage: map2.sh <clean_dir> <bam_dir> <hisat2|STAR>}
tool=${3:?usage: map2.sh <clean_dir> <bam_dir> <hisat2|STAR>}

HISAT2_INDEX=${HISAT2_INDEX:-/path/to/GRCm39/hisat2_index/Mus_musculus.GRCm39.dna.primary_assembly}
STAR_INDEX=${STAR_INDEX:-/path/to/GRCm39/star_index}
HISAT2_THREADS=${HISAT2_THREADS:-10}
STAR_THREADS=${STAR_THREADS:-30}

mkdir -p "$bam_dir"
shopt -s nullglob

case "$tool" in
hisat2)
    for fq1 in "$clean_dir"/*_1_val_1.fq.gz; do
        sample=$(basename "$fq1" _1_val_1.fq.gz)
        fq2="$clean_dir/${sample}_2_val_2.fq.gz"
        [ -f "$fq2" ] || { echo "ERROR: missing mate for ${sample}: $fq2" >&2; exit 1; }
        bam="$bam_dir/${sample}.Hisat_aln.sorted.bam"
        hisat2 -p "$HISAT2_THREADS" -x "$HISAT2_INDEX" -1 "$fq1" -2 "$fq2" \
            2> "$bam_dir/${sample}.hisat2.log" \
          | samtools sort -@ "$HISAT2_THREADS" -o "$bam" -
        samtools index "$bam"
    done
    ;;
STAR)
    for fq1 in "$clean_dir"/*_1_val_1.fq.gz; do
        sample=$(basename "$fq1" _1_val_1.fq.gz)
        fq2="$clean_dir/${sample}_2_val_2.fq.gz"
        [ -f "$fq2" ] || { echo "ERROR: missing mate for ${sample}: $fq2" >&2; exit 1; }
        STAR --runThreadN "$STAR_THREADS" \
             --genomeDir "$STAR_INDEX" \
             --readFilesCommand zcat \
             --readFilesIn "$fq1" "$fq2" \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "$bam_dir/${sample}.STAR."
        samtools index "$bam_dir/${sample}.STAR.Aligned.sortedByCoord.out.bam"
    done
    ;;
*)
    echo "ERROR: tool must be hisat2 or STAR" >&2
    exit 1
    ;;
esac

echo "Alignment done: $bam_dir"
