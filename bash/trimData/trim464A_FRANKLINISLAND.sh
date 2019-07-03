#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to) 
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 1            # number of cores
#SBATCH --time=6:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=4GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/trimData_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/trimData_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: script generate gene expression profiles for each venom sample

# Load  modules
module load AdapterRemoval

## Script for FastQC and adapter removal
ROOTDIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq/
DATADIR=${ROOTDIR}/0_rawData/fastq
F=${DATADIR}/464A_FRANKLINISLAND_R1.fastq.gz
OUT1=${ROOTDIR}/1_trimmedData/fastq/464A_FRANKLINISLAND.trimmed1.fq.gz

## Check the file exists
if [[ -f ${F} ]]; then
  echo "Found file ${F}"
fi

AdapterRemoval \
  --file1 ${F} \
  --file2 ${F%R1.fastq.gz}R2.fastq.gz \
  --gzip  \
  --trimns \
  --trimqualities \
  --minlength 20 \
  --output1 ${OUT1} \
  --output2 ${OUT1%1.fq.gz}2.fq.gz


