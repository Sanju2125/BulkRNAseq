
# Individual FASTQC Report
fastqc Rnaseq_Cancer/1.data/raw_fastq/* .fastq.gz \ 
-o Rnaseq_Cancer/2.qc_fastqc_raw/fastqc


# MULTIQC Report
multiqc Rnaseq_Cancer/2.qc_fastqc_raw \ 
-o Rnaseq_Cancer/2.qc_fastqc_raw/multiqc



