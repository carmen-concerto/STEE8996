module load uri/main
module load bwa/0.7.17
module load GATK/4.3.0.0-GCCcore-11.3.0-Java-11
module load R/4.3.2-gfbf-2023a
module load samtools/1.19.2
module load picard/3.1.1
module load freebayes/1.3.7-gfbf-2023a-R-4.3.2
module load bedtools2/2.31.1

A="A_CKDN220037208-1A_H3GJMDSX5_L3"

#bwa for mapping reads to reference
bwa mem -t 8 -M -a -v 1 $REF ${A}_1.fq.gz ${A}_2.fq.gz > bam/$A.sam

# Picard CleanSam
java -Xmx32g -jar $PICARD_JAR CleanSam I=bam/${A}.sam O=bam/${A}.clean.sam

#SAM to BAM
samtools view -b -T $REF -o bam/${A}.bam bam/${A}.clean.sam

#Picard FixMateInformation
java -Xmx32g -jar $PICARD_JAR FixMateInformation SO=coordinate I=bam/${A}.bam O=bam/${A}.fixmate.bam

#Picard MarkDuplicates
java -Xmx32g -jar $PICARD_JAR MarkDuplicates \
    I=bam/${A}.fixmate.bam \
    O=bam/${A}.fixmate.dedup.bam \
    M=dedup_metrics/${A}.fixmate.dedup.metrics

#Index BAM
samtools index bam/${A}.fixmate.dedup.bam

#Samtools Flagstat
samtools flagstat bam/${A}.fixmate.dedup.bam > flagstat/${A}.txt

module load Java/11.0.16
module load picard/2.25.1-Java-11

#Coverage stats
java -Xmx32g -jar "$PICARD_JAR" CollectRawWgsMetrics \
    I="bam/${A}.fixmate.dedup.bam" \
    O="wgs_metrics/${A}.raw_wgs_metrics" \
    R=$REF
    
java -Xmx10g -jar "$PICARD_JAR" CollectInsertSizeMetrics \
    I="bam/${A}.fixmate.dedup.bam" \
    O="insert_metrics/${A}.insertMetrics.txt" \
    H="insert_metrics/${A}.insertHistogram.pdf" \
    M=0.5 \
    ASSUME_SORTED=false

#FreeBayes variant calling
freebayes -f $REF -C 3 -p 1 "bam/${A}.fixmate.dedup.bam" > "fb_vcf/${A}.fb.vcf"

#GRIDSS structural variants
java -Xmx50g -cp "$GRIDSS_JAR" gridss.CallVariants \
    OUTPUT="GRIDSS/${A}.SV.vcf" \
    ASSEMBLY="GRIDSS/${A}.SV.assembly.bam" \
    INPUT="bam/${A}.fixmate.dedup.bam" \
    REFERENCE_SEQUENCE="$REF" \
    TMP_DIR="$TMP_DIR" \
    WORKING_DIR="$WORKING_DIR" \
    WORKER_THREADS=8 # Match SLURM ntasks
    
# Genome coverage
bedtools genomecov -d -ibam "bam/${A}.fixmate.dedup.bam" -g "${REF}.fai" > genomecov/${A}.genomecov.txt




