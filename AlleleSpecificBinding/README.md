# Context-TFs
Material and source code for Kribelbauer et al. (https://www.biorxiv.org/content/10.1101/2023.05.05.539543v1). 

# Allele Specific Binding scripts
For calculating Allele Specific Binding, we used a 3-step protocol:
1. Running bwa mem to align the fastq files to Homo_sapiens GRCh37.75 (hg19) assembly
2. Running Picard tool for removing duplicated reads
3. Running Freebayes with an input VCF file, to count reads mapping to each allele

We've run this for each TF separately. Designing different scripts depending on the protocol used: [ATAC-seq, ChIP-seq, GRO-seq, HyDRA] and [single-end, paired-end]

# Issues encountered
Freebayes has an inner bug, making it impossible to run the tool with both `--variant-input` and `--only-use-input-alleles` options activated. See GitHub issue (still unresolved to this day) here: (https://github.com/freebayes/freebayes/issues/166)[https://github.com/freebayes/freebayes/issues/166]
Therefore, we designed a work-around strategy by splitting the assembly by chromosome, and processing each chromosome separately.

