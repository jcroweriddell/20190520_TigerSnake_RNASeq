#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to)
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=6:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/featureCounts/fcounts_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/featureCounts/fcounts_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: Script to generate gene expression profiles for each venom sample

## Load modules
module load Subread/1.5.2-foss-2016b

## Set variables
ROOT=/fast/users/a1662801/20190520_TigerSnake_RNASeq
DATADIR=${ROOT}/2_alignedData/bam3
GTF=${ROOT}/4_assembly/stringtie_merged.gtf
OUTDIR=${ROOT}/5_readCounts

## Estimate read abundances using featureCounts
## Previous step: assemble transcripts using stringTie, merge output gtf files in single .gtf for all samples
# NOTE: -s 2 indicates forward reverse reads, use merged gtf from merged stringtie
featureCounts -T 8 \
	-s 2 \
	-g gene_name \
	-a ${GTF} \
	-o ${OUTDIR}/counts.txt \
	${DATADIR}/*out.bam

## Create matrix for downstream DE analyses
cut -f 1,7-52 ${OUTDIR}/counts.txt > ${OUTDIR}/countsMatrix.txt

# Clean up matrix for downstream DE analyses
# Delete first line of txt file & clean up sample names
sed -i '1d;$d' ${OUTDIR}/countsMatrix.txt
sed -i '1s:/fast/users/a1662801/20190520_TigerSnake_RNASeq/2_alignedData/bam3/:S:g' ${OUTDIR}/countsMatrix.txt
sed -i '1s:Aligned.sortedByCoord.out.bam::g' ${OUTDIR}/countsMatrix.txt


