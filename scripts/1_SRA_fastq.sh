set -euo pipefail

INPUT=scripts/srr_ids.txt
OUTDIR=data/raw_fastq
SRADIR=~/sra_cache/sra

mkdir -p "$OUTDIR"

while read -r sra; do
    echo "Converting $sra"
    fasterq-dump "$SRADIR/${sra}.sra" --split-files -O "$OUTDIR"
    gzip "$OUTDIR/${sra}_1.fastq"
    gzip "$OUTDIR/${sra}_2.fastq"
done < "$INPUT"
