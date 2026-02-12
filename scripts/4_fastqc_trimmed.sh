set -euo pipefail
#Individual FATSQC
fastqc data/trimmed_fastq/*.fastq.gz -o qc/fastqc_trimmed -t 8

#MULTIQC
multiqc qc/fastqc_trimmed -o qc/multiqc_trimmed
