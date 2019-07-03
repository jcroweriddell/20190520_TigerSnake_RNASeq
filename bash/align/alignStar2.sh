#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to)
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=48:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/align2_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/align2_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: Script to align new trimmed reads (min length set to 50 bp) using STAR to see if alignment rate is improved

# Load modules
#module load fastqc/0.11.4
#module load AdapterRemoval/2.2.1-foss-2016b
module load STAR/2.7.0d-foss-2016b
#module load StringTie/1.3.3-foss-2017a
#module load Subread/1.5.2-foss-2016b

# Set directory variables
ROOTDIR=/fast/users/a1662801/20190520_TigerSnake_RNASeq
#RAWREADS=${ROOTDIR}/0_rawData/fastq
TRIMREADS=${ROOTDIR}/1_trimmedData/fastq2
GENOME=/data/biorefs/reference_genomes/animals/Notechis_scutatus/star
BAMDIR=${ROOTDIR}/2_alignedData/bam4
#GTF=/fast/users/a1662801/20190520_TigerSnake_RNASeq/3_genomeRef/gtf/GCF_900518725.1_TS10Xv2-PRI_genomic.gtf
#ASSEMBLE=${ROOTDIR}/4_assembly
#COUNTDIR=${ROOTDIR}/5_readCounts


## Code for trimming and aligning RNA-seq reads
echo "Starting trimming and aligning reads"

for FQGZ in ${TRIMREADS}/*R1.trimmed1.fq.gz
    do
    # Find each read for paired end reads
    t1=${FQGZ}
    t2=${FQGZ%R1.trimmed1.fq.gz}R2.trimmed2.fq.gz
    echo "Read 1 is $(basename ${t1}) "
    echo "Read 2 is $(basename ${t2}) "

    # Align reads to reference genome using STAR
    # NOTE 2: need to have already build genome indices using STAR
    echo "Starting alignment of $(basename ${FQGZ%_R1.trimmed1.fq.gz}) "
    STAR --genomeDir ${GENOME} \
         --runThreadN 8 \
         --readFilesCommand gunzip -c \
         --readFilesIn ${t1} ${t2} \
         --outFileNamePrefix ${BAMDIR}/$(basename ${FQGZ%_R1.trimmed1.fq.gz}) \
         --outSAMtype BAM SortedByCoordinate
    echo "Finishing alignment of $(basename ${FQGZ%_R1.trimmed1.fq.gz}) "

done
