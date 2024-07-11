#!/bin/bash
#SBATCH --job-name="ASB_ENCODE3"
#SBATCH --error="ASB_ENCODE3_%A_%a.err"
#SBATCH --output="ASB_ENCODE3_%A_%a.out"
#SBATCH --nodes 1
#SBATCH --ntasks=1 # Total number of MPI tasks
#SBATCH --ntasks-per-node=1 # Number of MPI task per compute node (<=8)
#SBATCH --cpus-per-task=8 # 8
#SBATCH --time=12:00:00
#SBATCH --array=0-2506 # Job Array
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
TFs=(ASH2L ATF2 ATF3 BATF BCL11A BCL3 BCLAF1 BHLHE40 BRCA1 CBFB CBX3 CBX5 CEBPB CEBPZ CHD1 CHD2 CHD4 CREB1 CREM CTCF CUX1 E2F4 EBF1 EED EGR1 ELF1 ELK1 EP300 ESRRA ETS1 ETV6 EZH2 FOS FOXM1 GABPA HCFC1 HDAC2 HDAC6 HSF1 IKZF1 IRF3 IRF4 JUNB JUND KAT2A KDM1A MAFK MAX MAZ MEF2A MEF2C MTA3 MXI1 MYB MYC NFATC1 NFIC NFYA NFYB NR2C2 NRF1 PAX5 PBX3 PML POLR2A POLR3G POU2F2 RAD21 RBBP5 RCOR1 RELA REST RFX5 RUNX3 RXRA SIN3A SIX5 SMAD5 SMC3 SP1 SPI1 SREBF1 SREBF2 SRF STAT1 STAT3 STAT5A SUPT20H SUZ12 TAF1 TARDBP TBL1XR1 TBP TCF12 TCF3 TCF7 UBTF USF1 USF2 WRNIP1 YY1 ZBED1 ZBTB33 ZEB1 ZFP36 ZNF143 ZNF274 ZNF384 ZZZ3)
CHRs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
thresholds=(30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 20 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 20 30 30 30 30 20 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30)

# Current job in the job array
# 109 TFs x 23 CHRs = 2507 jobs: 0-2506
TF=${TFs[$((${SLURM_ARRAY_TASK_ID} / 23))]} # 000000000 and then 11111111111
thres=${thresholds[$((${SLURM_ARRAY_TASK_ID} / 23))]} # 000000000 and then 11111111111
chr=${CHRs[$((${SLURM_ARRAY_TASK_ID} % 23))]} # 0-23 and then 0-23 again

# Setting up paths
rootdir=/scratch/gardeux/ENCODE3_ASB_Judith
fastafile=${rootdir}/fasta/Homo_sapiens.GRCh37.75.dna.chromosome.${chr}.fa
vcffileIN=${rootdir}/input_vcf/GM12878_GRCh37_GIAB_highconf_noindels_chr${chr}.recode.vcf
bamfile=${rootdir}/bam/${TF}_chr${chr}
vcffileOUT=${rootdir}/regeno_vcf/${TF}_chr${chr}.regeno.vcf
fastqfile=${rootdir}/fastq/single-end/${TF}

# Create directory
mkdir -p ${rootdir}/bam/
mkdir -p ${rootdir}/regeno_vcf/

echo "Processing " ${TF} " with threshold " ${thres} " on chromosome " ${chr} " ...";

if [ -f "${fastqfile}.fastq.gz" ] 
then
	echo "1. Aligning with bwa";
	# -F 4 removes unmapped reads, -Sb transform to BAM, and sort
	# -T 20 only for small reads (default -T 30 is unreachable for reads <30nt)
	#${bwapath}/bwa mem -t 8 -T ${thres} -M -R '@RG\tID:'"$TF"'\tSM:'"$TF"'\tPL:illumina\tLB:'"$TF"'.library\tPU:'"$TF"'.unit' ${fastafile} ${fastqfile}.fastq.gz | samtools view -Sb -F 4 | samtools sort -o ${bamfile}.bam

	echo "2. Removing Duplicates";
	#java -jar ${picardpath}/picard.jar MarkDuplicates I=${bamfile}.bam O=${bamfile}.dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=${bamfile}.dedup.metrics.txt REMOVE_DUPLICATES=true

	echo "3. Regenotype all variants in BAM" # Removed --no-partial-observations
	${toolpath}/freebayes -f ${fastafile} --report-monomorphic --only-use-input-alleles --min-alternate-fraction 0 --variant-input ${vcffileIN} ${bamfile}.dedup.bam >${vcffileOUT}
else
    echo "No Control fastq files for " ${TF} ". Aborting..."
fi

exit 0;

