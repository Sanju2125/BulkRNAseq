set -euo pipefail

mkdir -p counts

featureCounts \
    -T 8 \
    -p \
    --countReadPairs \
    -B \
    -C \
    -t exon \
    -g gene_id \
    -a 4.reference/gencode.v44.annotation.gtf \
    -o counts/gene_counts.txt \
    alignment/bam/*.sorted.bam
