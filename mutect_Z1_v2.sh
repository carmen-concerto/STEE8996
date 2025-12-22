#!/bin/bash
#SBATCH --job-name=mutect2_Z1
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=mutect2_Z1_%j.out
#SBATCH --error=mutect2_Z1_%j.err

# Load modules
module load uri/main
module load Java/11.0.16
module load bwa/0.7.17
module load GATK/4.3.0.0-GCCcore-11.3.0-Java-11
module load samtools/1.19.2

# Reference genome
REF="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/ref_Dilay/Fusarium_oxysporum_MRL8996_dha.fasta"

# Tumor (evolved) BAM
TUMOR_BAM="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/Z1/A_CKDN220037208-1A_H3GJMDSX5_L3.fixmate.dedup.recal.bam"

# Normal (ancestral) BAM
NORMAL_BAM="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/D19/bam/D19_FDSW202340874-1r_H5HYHDSXY_L3.fixmate.dedup.RG.bam"
NORMAL_NAME="D19_FDSW202340874-1r_H5HYHDSXY_L3"

# Output files
OUT_VCF="Z1_vs_D19.vcf.gz"
FILTERED_VCF="Z1_vs_D19.filtered.vcf.gz"
F1R2="Z1.f1r2.tar.gz"

# Run Mutect2
gatk --java-options "-Xmx30g" Mutect2 \
  -R "$REF" \
  -I "$TUMOR_BAM" \
  -I "$NORMAL_BAM" \
  -normal "$NORMAL_NAME" \
  --tmp-dir "/scratch/workspace/siyuanwu_umass_edu-stee" \
  --genotype-germline-sites true \
  --f1r2-tar-gz "$F1R2" \
  -O "$OUT_VCF" \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --tumor-lod-to-emit 2.0 \
  --max-mnp-distance 0 \
  --pruning-lod-threshold 2.0

# Run FilterMutectCalls
gatk --java-options "-Xmx20g" FilterMutectCalls \
  -V "$OUT_VCF" \
  -R "$REF" \
  -O "$FILTERED_VCF"
