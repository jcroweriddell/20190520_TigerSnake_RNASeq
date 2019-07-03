#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to)
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=2:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/renameBAM_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/renameBAM_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: Renaming read alignment files

## Set variables
DATADIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq/2_alignedData/bam
OUTDIR=${DATADIR}/renameBamFiles

for B in ${DATADIR}/*out.bam
 do
    echo "moving $B to $OUTDIR and renaming to "$(basename ${B%.trimmed1Aligned.sortedByCoord.out.bam}).bam" "
    cp $B ${OUTDIR}/$(basename ${B%.trimmed1Aligned.sortedByCoord.out.bam}).bam 
    echo "finished with $B"
done
