#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to)
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=16:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/alignStar/alignStar_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/alignStar/alignStar_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: Script to generate gene expression profiles for each venom sample

# Load modules
module load STAR/2.7.0d-foss-2016b

# Set variables
DATADIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq/1_trimmedData/fastq
OUTDIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq/2_alignedData/bam2

## Align reads using STAR
 
## Previous step: build genome indices using STAR generate genome and annotation gtf file

for FQGZ in ${DATADIR}/*trimmed1.fq.gz
 do
    # Find each read for paired end reads
    f1=${FQGZ}
    f2=${FQGZ%trimmed1.fq.gz}trimmed2.fq.gz
    echo "Read 1 is $f1"
    echo "Read 2 is $f2"

    # Align trimmed reads against reference using STAR
    echo "Starting alignment of "$(basename ${FQGZ%.trimmed1.fq.gz})" "
    
    STAR --genomeDir /data/biorefs/reference_genomes/animals/Notechis_scutatus/star \
         --runThreadN 8 \
         --readFilesCommand gunzip -c \
         --readFilesIn ${f1} ${f2} \
         --outFileNamePrefix ${OUTDIR}/$(basename ${FQGZ%.trimmed1.fq.gz}) \
         --outSAMtype BAM SortedByCoordinate 
    echo "Finishing alignment of "$(basename ${FQGZ%.trimmed1.fq.gz})" "
done

