#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to) 
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=8:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/buildIndex_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/buildIndex_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: script generate gene expression profiles for each venom sample


# Load our modules
module load HISAT2/2.0.5-foss-2016uofa
module load Java/1.8.0_71


## Build reference index
# Using hisat2
<<<<<<< HEAD
hisat2-build /fast/users/a1662801/20190520_TigerSnake_RNASeq/3_genomeRef/fastq/GCF_900518725.1_TS10Xv2-PRI_genomic.fna tigerSnakeRef
=======
hisat2-build /fast/users/a1662801/20190520_TigerSnake_RNASeq/3_genomeRef/fastq/GCF_900518725.1_TS10Xv2-PRI_genomic.fna.gz tigerSnakeRef
>>>>>>> 5ec0f217af2dd1ae3ee30afab9f33136dd63ff5c

