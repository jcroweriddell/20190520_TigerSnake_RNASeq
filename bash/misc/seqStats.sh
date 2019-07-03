#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to)
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=6:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/stats_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/stats_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: Script to generate gene expression profiles for each venom sample

## Load modules
module load seqkit/0.8.1

## Set variables
ROOT=/fast/users/a1662801/20190520_TigerSnake_RNASeq
OUTDIR=${ROOT}/stats

# Run sequencing statics on raw reads and trimmed reads
seqkit stats ${ROOT}/0_rawData/fastq/*.fastq.gz > ${OUTDIR}/raw2.stats
seqkit stats ${ROOT}/1_trimmedData/fastq/*.fq.gz > ${OUTDIR}/trim2.stats
