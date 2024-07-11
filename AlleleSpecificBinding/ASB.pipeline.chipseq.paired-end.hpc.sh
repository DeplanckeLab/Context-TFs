#!/bin/bash
#SBATCH --job-name="ASB_ENCODE3"
#SBATCH --error="ASB_ENCODE3_%A_%a.err"
#SBATCH --output="ASB_ENCODE3_%A_%a.out"
#SBATCH --nodes 1
#SBATCH --ntasks=1 # Total number of MPI tasks
#SBATCH --ntasks-per-node=1 # Number of MPI task per compute node (<=8)
#SBATCH --cpus-per-task=8 # 8
#SBATCH --time=12:00:00
#SBATCH --array=0-1356 # Job Array
#SBATCH --mem=32G # per node. Or use --mem-per-cpu 
#SBATCH --mail-type=NONE
#SBATCH --mail-user=vincent.gardeux@epfl.ch

# Load required modules from scitas
module load intel/18.0.5 samtools/1.9 python/3.7.3 fastqc/0.11.8

# Load specific tools that we installed ourselves
toolpath=/home/gardeux/Tools
deeptoolspath=${toolpath}/deepTools-3.2.1/bin
picardpath=${toolpath}/picard-2.17.8
bwapath=${toolpath}/bwa-0.7.17
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${toolpath}:${deeptoolspath}:${bwapath}:${picardpath}
PATH=${PATH}:${toolpath}:${deeptoolspath}:${bwapath}:${picardpath}

# All jobs infos
TFs=(ARID3A ARNT ATF2 ATF7 BACH1 BCLAF1 BHLHE40 BMI1 DPF2 E2F8 E4F1 ELF1 ETV6 FOXK2 GATAD2B HDGF IKZF1 IKZF2 IRF3 IRF5 KLF5 LARP7 MAZ MEF2B MLLT1 MTA2 NBN NFATC3 NFXL1 NKRF NR2C1 NR2F1 PAX8 PKNOX1 PRDM15 RAD51 RB1 RELB SKIL SMAD1 SMARCA5 SRF STAT1 TARDBP TBX21 TCF12 TRIM22 YBX1 ZBTB17 ZBTB33 ZBTB40 ZNF143 ZNF207 ZNF217 ZNF24 ZNF592 ZNF622 ZNF687 ZSCAN29)
CHRs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
thresholds=(30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30)

# Current job in the job array
# 59 TFs x 23 CHRs = 1357 jobs: 0-1356
TF=${TFs[$((${SLURM_ARRAY_TASK_ID} / 23))]} # 000000000 and then 11111111111
thres=${thresholds[$((${SLURM_ARRAY_TASK_ID} / 23))]} # 000000000 and then 11111111111
chr=${CHRs[$((${SLURM_ARRAY_TASK_ID} % 23))]} # 0-23 and then 0-23 again

# Setting up paths
rootdir=/scratch/gardeux/ENCODE3_ASB_Judith
fastafile=${rootdir}/fasta/Homo_sapiens.GRCh37.75.dna.chromosome.${chr}.fa
vcffileIN=${rootdir}/input_vcf/GM12878_GRCh37_GIAB_highconf_noindels_chr${chr}.recode.vcf
bamfile=${rootdir}/bam/${TF}_PE_chr${chr}
vcffileOUT=${rootdir}/regeno_vcf/${TF}_PE_chr${chr}.regeno.vcf
fastqfile=${rootdir}/fastq/paired-end/${TF}

# Create directory
mkdir -p ${rootdir}/bam/
mkdir -p ${rootdir}/regeno_vcf/

echo "Processing " ${TF} " with threshold " ${thres} " on chromosome " ${chr} " ...";

if [ -f "${fastqfile}_R1.fastq.gz" ] 
then
	echo "1. Aligning with bwa";
	# -F 4 removes unmapped reads, -Sb transform to BAM, and sort
	# -T 20 only for small reads (default -T 30 is unreachable for reads <30nt)
	${bwapath}/bwa mem -t 8 -T ${thres} -M -R '@RG\tID:'"$TF"'\tSM:'"$TF"'\tPL:illumina\tLB:'"$TF"'.library\tPU:'"$TF"'.unit' ${fastafile} ${fastqfile}_R1.fastq.gz ${fastqfile}_R2.fastq.gz | samtools view -Sb -F 4 | samtools sort -o ${bamfile}.bam

	echo "2. Removing Duplicates";
	java -jar ${picardpath}/picard.jar MarkDuplicates I=${bamfile}.bam O=${bamfile}.dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=${bamfile}.dedup.metrics.txt REMOVE_DUPLICATES=true

	echo "3. Regenotype all variants in BAM" # Removed --no-partial-observations
	${toolpath}/freebayes -f ${fastafile} --report-monomorphic --only-use-input-alleles --min-alternate-fraction 0 --variant-input ${vcffileIN} ${bamfile}.dedup.bam >${vcffileOUT}
else
    echo "No Control R1 fastq files for " ${TF} ". Aborting..."
fi

exit 0;
