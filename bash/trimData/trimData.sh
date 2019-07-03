#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to) 
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=1:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/trimData_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/trimData_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: script generate gene expression profiles for each venom sample

# Load  modules
module load AdapterRemoval

## Script for FastQC and adapter removal
data_dir=/fast/users/a1662801/20190520_TigerSnake_RNASeq/0_rawData/fastqRemaining

AdapterRemoval --file1 $data_dir/R14_REEVESBYISLAND_R1.fastq.gz --file2 $data_dir/R14_REEVESBYISLAND_R2.fastq.gz \
--gzip  --trimns --trimqualities --minlength 20

