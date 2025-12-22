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


TEFINDER_PATH="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/scripts/TEfinder"

REF="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/ref_Dilay/Fusarium_oxysporum_MRL8996_dha.fasta"
TE_gff="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/TEfinder/Fusarium_oxysporum_MRL8996_JGI_transposons.gff"
TE_list="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/TEfinder/TE_list_clean.txt"


D19_bam="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/D19/bam/D19_FDSW202340874-1r_H5HYHDSXY_L3.fixmate.dedup.RG.bam"
Z1_bam="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/Z1/bam/A_CKDN220037208-1A_H3GJMDSX5_L3.fixmate.dedup.RG.bam"
Z2_bam="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/Z2/bam/B1_CKDN220039295-1A_H37TKDSX5_L2.fixmate.dedup.RG.bam"
Z3_bam="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/Z3/bam/C_CKDN220037210-1A_H3GJMDSX5_L3.fixmate.dedup.RG.bam"
Z4_bam="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/Z4/bam/D_CKDN220037211-1A_H3GJMDSX5_L3.fixmate.dedup.RG.bam"
Z5_bam="/work/pi_lijun_umass_edu/users/Melanie/projects/STEE_mouse/STEE_mouse/Z5/bam/E_CKDN220037212-1A_H3GJMDSX5_L3.fixmate.dedup.RG.bam"

#mkdir -p TEfinder_out_D19

date
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

echo "Finished Z1"
date

#Z2
#${TEFINDER_PATH}/TEfinder \
#  -alignment ${Z2_bam} \
#  -fa ${REF} \
#  -gtf ${TE_gff} \
#  -te ${TE_list} \
#  -fis 500 \
#  -k 30 \
#  -out GTF \
#  -intermed yes \
#  --pathToPicardjar ${PICARD_JAR} \
#  -workingdir TEfinder_out_Z2_${SLURM_JOB_ID} \
#  -outname Z2
#
#echo "Finished Z2"
#date

#Z3
#${TEFINDER_PATH}/TEfinder \
#  -alignment ${Z3_bam} \
#  -fa ${REF} \
#  -gtf ${TE_gff} \
#  -te ${TE_list} \
#  -fis 500 \
#  -k 30 \
#  -out GTF \
#  -intermed yes \
#  --pathToPicardjar ${PICARD_JAR} \
#  -workingdir TEfinder_out_Z3_${SLURM_JOB_ID} \
#  -outname Z3

#echo "Finished Z3"
#date

#Z4
#${TEFINDER_PATH}/TEfinder \
#  -alignment ${Z4_bam} \
#  -fa ${REF} \
#  -gtf ${TE_gff} \
#  -te ${TE_list} \
#  -fis 500 \
#  -k 30 \
#  -out GTF \
#  -intermed yes \
#  --pathToPicardjar ${PICARD_JAR} \
#  -workingdir TEfinder_out_Z4_${SLURM_JOB_ID} \
#  -outname Z4

#echo "Finished Z4"
#date

#Z5
#${TEFINDER_PATH}/TEfinder \
#  -alignment ${Z5_bam} \
#  -fa ${REF} \
#  -gtf ${TE_gff} \
#  -te ${TE_list} \
#  -fis 500 \
#  -k 30 \
#  -out GTF \
#  -intermed yes \
#  --pathToPicardjar ${PICARD_JAR} \
#  -workingdir TEfinder_out_Z5_${SLURM_JOB_ID} \
#  -outname Z5

#echo "Finished Z5"
#date
