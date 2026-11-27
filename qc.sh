#!/bin/bash
# Read-level QC of raw FASTQ files: FastQC per file, then MultiQC over the set.
#
#   bash qc.sh <fastq_dir> [qc_dir]
#
# Input   <fastq_dir>/*.fq.gz or *.fastq.gz
# Output  <qc_dir> (default <fastq_dir>/qc): FastQC html/zip plus multiqc_report.html
#
# fastqc and multiqc must be on PATH; no conda environment is activated here.
set -euo pipefail

fq_dir=${1:?usage: qc.sh <fastq_dir> [qc_dir]}
qc_dir=${2:-${fq_dir%/}/qc}

mkdir -p "$qc_dir"
shopt -s nullglob
files=("$fq_dir"/*.fq.gz "$fq_dir"/*.fastq.gz)
[ ${#files[@]} -gt 0 ] || { echo "ERROR: no .fq.gz / .fastq.gz in $fq_dir" >&2; exit 1; }

echo "FastQC on ${#files[@]} file(s)"
fastqc -t 8 -o "$qc_dir" "${files[@]}"
multiqc "$qc_dir" -o "$qc_dir"

echo "QC reports: $qc_dir"
