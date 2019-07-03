#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to)
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=6:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=16GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/stringTie/stringtie2_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/stringTie/stringtie2_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: Script to generate gene expression profiles for each venom sample

# Load modules
module load StringTie/1.3.3-foss-2017a

# Set directory variables
DATADIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq/2_alignedData/bam3
GTF=/fast/users/a1662801/20190520_TigerSnake_RNASeq/3_genomeRef/gtf/GCF_900518725.1_TS10Xv2-PRI_genomic.gtf
OUTDIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq/4_assembly

## Assemble RNA-seq alignments into transcripts using Stringtie
## Previous step: align trimmed reads to reference genome using STAR

# set count variable to count total number of files that have been assembled in loop
count=0

for B in ${DATADIR}/*out.bam
 do
    # Assemble transcripts
    echo "Start assembling transcripts for "$(basename ${B%Aligned.sortedByCoord.out.bam})" "
    stringtie $B \
        -l $(basename ${B%Aligned.sortedByCoord.out.bam}) \
        -p 8 \
        -G ${GTF} \
        -o ${OUTDIR}/"$(basename ${B%Aligned.sortedByCoord.out.bam})".gtf \
        -v 
    echo "Finish assembling transcripts for "$(basename ${B%Aligned.sortedByCoord.out.bam})" "
    ((count++))
    echo " $count out of 42 complete "
done

## Merge transcripts into single gtf file
# Create a mergelist.txt file that has list of each gtf file for each sample
echo "Create a mergelist from stringTie gtf output files"

find ${OUTDIR} -mindepth 1 -type f -name "*.gtf" > ${OUTDIR}/mergelist.txt

echo "Merge all gtf files"

# Merge all transcripts
stringtie --merge \
    -p 8 \
    -G ${GTF} \
    -o ${OUTDIR}/stringtie_merged.gtf ${OUTDIR}/mergelist.txt

echo "Merging of gtf files complete"

