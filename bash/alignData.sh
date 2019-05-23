#!/bin/bash
#SBATCH -p batch        # partition (this is the queue your job will be added to) 
#SBATCH -N 1            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8            # number of cores
#SBATCH --time=24:00:00 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=32GB      # memory pool for all cores (here set to 32 GB)
#SBATCH -o /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/trimData_%j.out
#SBATCH -e /fast/users/a1662801/20190520_TigerSnake_RNASeq/slurm/trimData_%j.err
#SBATCH --mail-type=END    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=jenna.crowe-riddell@adelaide.edu.au  # Email to which notification will be sent

## DESCRIPTION: script generate gene expression profiles for each venom sample


# Load our modules
module load HISAT2/2.0.5-foss-2016uofa
module load HTSeq/0.6.1p1-intel-2015c-Python-2.7.11
module load Java/1.8.0_71


## Build reference index
# Using hisat2
hisat2-build /fast/users/a1662801/20190520_TigerSnake_RNASeq/3_genomeRef/fastq/GCF_900518725.1_TS10Xv2-PRI_genomic.fna.gz
 tigerSnakeGenome

## Script for aligning trimmed reads to tiger snake reference genome
working=$(pwd)

data_dir=/fast/users/a1662801/20190520_TigerSnake_RNASeq/1_trimmedData/fastq
output_dir=/20190520_TigerSnake_RNASeq/2_alignedData/
genome_prefix=/fast/users/a1662801/20190520_TigerSnake_RNASeq/3_genomeRef/fastq/GCF_900518725.1_TS10Xv2-PRI_genomic.fna.gz

# Tiger snake GFF annotation file
gff=/fast/users/a1662801/20190520_TigerSnake_RNASeq/3_genomeRef/gff/GCF_900518725.1_TS10Xv2-PRI_genomic.gff.gz

# mkdir -p $output_dir

# Run fastqc on each read pairs
echo "starting fastqc"

fastqc -o $data_dir/FastQC $data_dir/fastq/*.fastq.gz 

echo "finished fastqc"

for FQGZ in $data_dir/fastq/*R1*.fastq.gz
 do
    # Everything indented (and before the "done") will be run
    
    # Align our trimmed reads against our reference
    echo "Starting Alignment of "$FQGZ" "
    hisat2 -x $genome_prefix \
	-1 $output_dir/fastq/"$(basename $FQGZ _R1.fastq.gz)".trimmed1.fq.gz \
	-2 $output_dir/fastq/"$(basename $FQGZ _R1.fastq.gz)".trimmed2.fq.gz | \
	samtools view -bS - > $output_dir/"$(basename $FQGZ _R1.fastq.gz)".hisat2Align.bam
    echo "Finishing Alignment of "$FQGZ" "

    # Stats
    samtools flagstat $output_dir/"$(basename $FQGZ _R1.fastq.gz)".hisat2Align.bam
    
    # sort before HTSeq	
    samtools sort -n $output_dir/"$(basename $FQGZ _R1.fastq.gz)".hisat2Align.bam \
	$output_dir/"$(basename $FQGZ _R1.fastq.gz)".hisat2Align.sorted

    # HTSeq (THIS IS NOT FPKM!!! Its Raw counts that are NOT normalised)
    htseq-count -f bam  "$(basename $FQGZ _R1.fastq.gz)".hisat2Algin.sorted.bam $gff
   
done
