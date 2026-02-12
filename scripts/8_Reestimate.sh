set -euo pipefail

SAMPLES=$(cat scripts/srr_ids.txt)

for s in $SAMPLES; do
    stringtie alignment/bam/${s}.sorted.bam \
        -p 8 \
        -G 6.stringtie/stringtie_merged.gtf \
        -e -B \
        -o 6.stringtie/${s}.merged.gtf
done
