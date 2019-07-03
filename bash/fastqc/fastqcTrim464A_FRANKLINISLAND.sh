#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to) 
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 1            # number of cores
#SBATCH --time=3:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=4GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/fastqcTrim_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/fastqcTrim_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: script generate gene expression profiles for each venom sample

# Set variables
DATADIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq/1_trimmedData/fastq
OUTDIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq/1_trimmedData/FastQC

mkdir ${OUTDIR}/fastq464A_FRANKLINISLAND

# Load modules
module load fastqc/0.11.4

# Run fastqc on each read pairs
echo "starting fastqc"

fastqc -o ${OUTDIR}/fastq464A_FRANKLINISLAND ${DATADIR}/464A_FRANKLINISLAND.trimmed* 

echo "finished fastqc"


