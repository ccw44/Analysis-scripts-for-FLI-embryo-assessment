#!/bin/bash
# Index the BAM files written by map2.sh, whichever aligner was used.
#
#   bash bai.sh <bam_dir>
set -euo pipefail

bam_dir=${1:?usage: bai.sh <bam_dir>}
shopt -s nullglob

count=0
for bam in "$bam_dir"/*.Hisat_aln.sorted.bam \
           "$bam_dir"/*.STAR.Aligned.sortedByCoord.out.bam; do
    if [ -f "${bam}.bai" ]; then
        echo "exists: ${bam}.bai"
        continue
    fi
    samtools index "$bam"
    echo "indexed: $bam"
    count=$((count + 1))
done

echo "${count} BAM file(s) indexed"
