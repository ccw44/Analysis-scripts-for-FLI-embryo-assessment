#!/bin/bash
# Adapter and quality trimming of paired-end reads with Trim Galore.
#
#   bash cleanfq.sh <fastq_dir> <project_dir> [cores_per_pair]
#
# Input   <fastq_dir>/<sample>_1.fq.gz and <fastq_dir>/<sample>_2.fq.gz
# Output  <project_dir>/cleanfq/<sample>_1_val_1.fq.gz
#         <project_dir>/cleanfq/<sample>_2_val_2.fq.gz
#
# The two mate names are enforced, so a sample missing its pair stops the run
# instead of producing a half-processed set. The trimmed file names are the ones
# map2.sh expects.
set -euo pipefail

fq_dir=${1:?usage: cleanfq.sh <fastq_dir> <project_dir> [cores_per_pair]}
proj_dir=${2:?usage: cleanfq.sh <fastq_dir> <project_dir> [cores_per_pair]}
cores=${3:-2}

clean_dir="${proj_dir%/}/cleanfq"
mkdir -p "$clean_dir"
shopt -s nullglob

pairs=0
for fq1 in "$fq_dir"/*_1.fq.gz; do
    sample=$(basename "$fq1" _1.fq.gz)
    fq2="$fq_dir/${sample}_2.fq.gz"
    if [ ! -f "$fq2" ]; then
        echo "ERROR: missing mate for ${sample}: $fq2" >&2
        exit 1
    fi
    # parameters as used in the original run
    trim_galore -q 25 --phred33 --length 36 --stringency 3 --paired \
                --cores "$cores" -o "$clean_dir" "$fq1" "$fq2"
    pairs=$((pairs + 1))
done

[ "$pairs" -gt 0 ] || { echo "ERROR: no *_1.fq.gz in $fq_dir" >&2; exit 1; }
echo "Trimmed ${pairs} pair(s) into $clean_dir"
