#!/bin/bash
#SBATCH --job-name=A_preGATK
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8   # Adjust based on resources
#SBATCH --mem=55G    # Adjust memory allocation
#SBATCH --time=12:00:00
#SBATCH --output=A_preGATK_%j.out
#SBATCH --error=A_preGATK_%j.err

#load modules

module load uri/main
module load bwa/0.7.17
module load GATK/4.3.0.0-GCCcore-11.3.0-Java-11
module load R/4.3.2-gfbf-2023a
module load samtools/1.19.2
module load picard/3.1.1
module load freebayes/1.3.7-gfbf-2023a-R-4.3.2
module load bedtools2/2.31.1

A="A_CKDN220037208-1A_H3GJMDSX5_L3"
REF="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/ref_Dilay/Fusarium_oxysporum_MRL8996_dha.fasta"
PICARD_JAR="/modules/spack/packages/linux-ubuntu24.04-x86_64_v3/gcc-13.2.0/picard-3.1.1-bt7ub4b3bxdokkur7tgujfx45oq7zjcj/bin/picard.jar"
GRIDSS_JAR="/work/pi_lijun_umass_edu/users/Melanie/software/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar"
TMP_DIR="/scratch3/workspace/siyuanwu_umass_edu-bqsr_temp"
WORKING_DIR="/scratch3/workspace/siyuanwu_umass_edu-bqsr_temp"

#indexing the reference files
#date
#echo “Indexing the reference file”
#samtools faidx $REF
#gatk CreateSequenceDictionary -R $REF -O "/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/ref_8996/MRL8996_ref.dict"
#bwa index $REF

date

# Create output directories if they don't exist
mkdir -p bam dedup_metrics flagstat

#bwa for mapping reads to reference
echo “bwa for mapping reads to reference”
bwa mem -t 8 -M -a -v 1 $REF ${A}_1.fq.gz ${A}_2.fq.gz > bam/$A.sam
echo

#Run Picard CleanSam

date
echo "Running Picard CleanSam..."
java -Xmx32g -jar $PICARD_JAR CleanSam I=bam/${A}.sam O=bam/${A}.clean.sam
date

# Convert SAM to BAM
date
echo "Running Samtools View..."
samtools view -b -T $REF -o bam/${A}.bam bam/${A}.clean.sam

date
echo "Samtools conversion complete."

# Run Picard FixMateInformation
date
echo "Running Picard FixMateInformation..."
java -Xmx32g -jar $PICARD_JAR FixMateInformation SO=coordinate I=bam/${A}.bam O=bam/${A}.fixmate.bam

# Run Picard MarkDuplicates
date
echo "Running Picard MarkDuplicates..."
java -Xmx32g -jar $PICARD_JAR MarkDuplicates \
    I=bam/${A}.fixmate.bam \
    O=bam/${A}.fixmate.dedup.bam \
    M=dedup_metrics/${A}.fixmate.dedup.metrics

date
echo "Picard MarkDuplicates complete."

# Index BAM using Samtools
date
echo "Running Samtools Index..."
samtools index bam/${A}.fixmate.dedup.bam

# Run Samtools Flagstat
date
echo "Running Samtools Flagstat..."
samtools flagstat bam/${A}.fixmate.dedup.bam > flagstat/${A}.txt

date

# loading more modules and downgrade picard for the following tasks
module load Java/11.0.16
module load picard/2.25.1-Java-11

PICARD_JAR="$EBROOTPICARD/picard.jar"

mkdir -p wgs_metrics

# Get coverage statistics
date
echo "Running Picard CollectRawWgsMetrics..."
java -Xmx32g -jar "$PICARD_JAR" CollectRawWgsMetrics \
    I="bam/${A}.fixmate.dedup.bam" \
    O="wgs_metrics/${A}.raw_wgs_metrics" \
    R=$REF

date

# Create output directories if they don't exist
mkdir -p insert_metrics fb_vcf SV

# Plot insert size histogram
date
echo "Running Picard CollectInsertSizeMetrics..."
java -Xmx10g -jar "$PICARD_JAR" CollectInsertSizeMetrics \
    I="bam/${A}.fixmate.dedup.bam" \
    O="insert_metrics/${A}.insertMetrics.txt" \
    H="insert_metrics/${A}.insertHistogram.pdf" \
    M=0.5 \
    ASSUME_SORTED=false


# Run FreeBayes for variant calling
date
echo "Running FreeBayes..."
freebayes -f $REF -C 3 -p 1 "bam/${A}.fixmate.dedup.bam" > "fb_vcf/${A}.fb.vcf"

date

# Create output directories if they don't exist
mkdir -p GRIDSS

# Run GRIDSS
date
echo "Running GRIDSS..."
java -Xmx50g -cp "$GRIDSS_JAR" gridss.CallVariants \
    OUTPUT="GRIDSS/${A}.SV.vcf" \
    ASSEMBLY="GRIDSS/${A}.SV.assembly.bam" \
    INPUT="bam/${A}.fixmate.dedup.bam" \
    REFERENCE_SEQUENCE="$REF" \
    TMP_DIR="$TMP_DIR" \
    WORKING_DIR="$WORKING_DIR" \
    WORKER_THREADS=8 # Match SLURM ntasks
date

mkdir -p genomecov

# Generate genome coverage
date
echo "Running bedtools genomecov..."
bedtools genomecov -d -ibam "bam/${A}.fixmate.dedup.bam" -g "${REF}.fai" > genomecov/${A}.genomecov.txt
date

echo “Pipeline finished running”
