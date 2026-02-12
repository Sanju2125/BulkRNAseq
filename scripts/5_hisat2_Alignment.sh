set -euo pipefail

SAMPLES=$(cat scripts/srr_ids.txt)

for s in $SAMPLES; do
    echo "Aligning $s"

    hisat2 -p 8 \
        --dta \
        -x 4.reference/GRCh38_hisat2 \
        -1 data/trimmed_fastq/${s}_1.trimmed.fastq.gz \
        -2 data/trimmed_fastq/${s}_2.trimmed.fastq.gz \
        2> alignment/stats/${s}.hisat2.log \
    | samtools sort -@ 8 -o alignment/bam/${s}.sorted.bam

done
