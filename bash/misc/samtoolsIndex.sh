#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to)
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=6:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=16GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/samtoolsIndex_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/samtoolsIndex_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: Script to generate gene expression profiles for each venom sample

# Load modules
module load SAMtools/1.9-foss-2016b

# Set variables
DATADIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq/2_alignedData/bam3

# Create a index .bai files from .bam alignment files

for B in ${DATADIR}/*out.bam
	do echo "Creating index file for $B"
	samtools index $B 
	echo "Finishing index file for $B"
done
