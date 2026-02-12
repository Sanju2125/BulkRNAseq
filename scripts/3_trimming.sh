set -euo pipefail

SAMPLES=$(cat scripts/srr_ids.txt)

mkdir -p data/trimmed_fastq
mkdir -p qc/fastp_reports

for s in $SAMPLES; do
    fastp \
        -i data/raw_fastq/${s}_1.fastq.gz \
        -I data/raw_fastq/${s}_2.fastq.gz \
        -o data/trimmed_fastq/${s}_1.trimmed.fastq.gz \
        -O data/trimmed_fastq/${s}_2.trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --cut_tail \
        --cut_mean_quality 20 \
        --length_required 30 \
        --thread 8 \
        --html qc/fastp_reports/${s}.html
done
