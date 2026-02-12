set -euo pipefail

SAMPLES=$(cat scripts/srr_ids.txt)

mkdir -p 6.stringtie

for s in $SAMPLES; do
    stringtie alignment/bam/${s}.sorted.bam \
        -p 8 \
        -G 4.reference/gencode.v44.annotation.gtf \
        -o 6.stringtie/${s}.gtf \
        -l ${s}
done
