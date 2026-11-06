#!/bin/bash
#SBATCH --job-name=tefinder_Z1
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8   # Adjust based on resources
#SBATCH --mem=80G    # Adjust memory allocation
#SBATCH --time=02:00:00
#SBATCH --output=tefinder_Z1_%j.out
#SBATCH --error=tefinder_Z1_%j.err


module load samtools/1.19.2
module load bedtools2/2.31.1
module load picard/3.1.1

PICARD_JAR="/modules/spack/packages/linux-ubuntu24.04-x86_64_v3/gcc-13.2.0/picard-3.1.1-bt7ub4b3bxdokkur7tgujfx45oq7zjcj/bin/picard.jar"

TEFINDER_PATH="TEfinder"

REF="Fusarium_oxysporum_MRL8996_dha.fasta"
TE_gff="Fusarium_oxysporum_MRL8996_JGI_transposons.gff"
TE_list="TE_list_clean.txt"


D19_bam="D19_FDSW202340874-1r_H5HYHDSXY_L3.fixmate.dedup.RG.bam"
Z1_bam="A_CKDN220037208-1A_H3GJMDSX5_L3.fixmate.dedup.RG.bam"

#Z1
${TEFINDER_PATH}/TEfinder \
  -alignment ${Z1_bam} \
  -fa ${REF} \
  -gtf ${TE_gff} \
  -te ${TE_list} \
  -fis 500 \
  -k 30 \
  -out GTF \
  -intermed yes \
  -maxHeapMem 16000 \
  --pathToPicardjar ${PICARD_JAR} \
  -workingdir TEfinder_out_Z1_${SLURM_JOB_ID} \
  -outname Z1
