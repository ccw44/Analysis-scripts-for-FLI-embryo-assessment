#!/bin/bash
# featureCounts gene-level and transcript-level counts on the BAMs from map2.sh.
#
#   bash count.sh <bam_dir> <out_dir> <hisat2|STAR>
#
# Output  <out_dir>/gene_counts_raw.txt        full featureCounts output
#         <out_dir>/Gene_Counts_Clean.csv      Geneid, Chr, samples, gene_name
#         <out_dir>/transcript_counts_raw.txt
#         <out_dir>/transcript_counts_clean.csv
#         <out_dir>/sampleid.txt
#
# GTF must be the same release used to build the index in map2.sh.
#
# Note on -p: it counts fragments. subread >= 2.0.2 additionally requires
# --countReadPairs for that meaning; check the subread version before comparing
# with the original counts (see RNAseq/README.md).
set -euo pipefail

bam_dir=${1:?usage: count.sh <bam_dir> <out_dir> <hisat2|STAR>}
out_dir=${2:?usage: count.sh <bam_dir> <out_dir> <hisat2|STAR>}
tool=${3:?usage: count.sh <bam_dir> <out_dir> <hisat2|STAR>}

GTF=${GTF:-/path/to/GRCm39/Mus_musculus.GRCm39.111.gtf}
THREADS=${THREADS:-20}

mkdir -p "$out_dir"
shopt -s nullglob

case "$tool" in
    hisat2) bam_suffix=".Hisat_aln.sorted.bam" ;;
    STAR)   bam_suffix=".STAR.Aligned.sortedByCoord.out.bam" ;;
    *)      echo "ERROR: tool must be hisat2 or STAR" >&2; exit 1 ;;
esac

: > "$out_dir/sampleid.txt"
for bam in "$bam_dir"/*"$bam_suffix"; do
    name=$(basename "$bam")
    echo "${name%"$bam_suffix"}" >> "$out_dir/sampleid.txt"
done
echo "samples: $(wc -l < "$out_dir/sampleid.txt")"

featureCounts -T "$THREADS" -p \
    -a "$GTF" -t exon -g gene_id --extraAttributes gene_name \
    -o "$out_dir/gene_counts_raw.txt" \
    "$bam_dir"/*"$bam_suffix"

sed '1d' "$out_dir/gene_counts_raw.txt" | cut -f1,2,7- | tr '\t' ',' \
    > "$out_dir/Gene_Counts_Clean.csv"

featureCounts -T "$THREADS" -p \
    -a "$GTF" -t exon -g transcript_id --extraAttributes gene_id,gene_name \
    -o "$out_dir/transcript_counts_raw.txt" \
    "$bam_dir"/*"$bam_suffix"

sed '1d' "$out_dir/transcript_counts_raw.txt" | cut -f1,2,7- | tr '\t' ',' \
    > "$out_dir/transcript_counts_clean.csv"

echo "Counts written to $out_dir"
