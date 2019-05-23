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


# Load  modules
module load fastqc/0.11.4
module load AdapterRemoval/2.2.0-foss-2016uofa

## Script for FastQC and adapter removal
working=$(pwd)

data_dir=/fast/users/a1662801/20190520_TigerSnake_RNASeq/0_rawData/
output_dir=/20190520_TigerSnake_RNASeq/1_trimmedData/

# Run FastQC on each read pair
echo "Starting fastqc"

fastqc -o $data_dir/FastQC $data_dir/fastq/*.fastq.gz 

echo "Finished fastqc"

for FQGZ in $data_dir/fastq/*R1*.fastq.gz
 do
    # Everything indented (and before the "done") will be run
 
    # Get our raw data and trim
    echo "Starting trimming of "$FQGZ" "
    AdapterRemoval --file1 $FQGZ --file2 ${FQGZ/R1/R2} \
	 --output1  $output_dir/fastq/"$(basename $FQGZ _R1.fastq.gz)".trimmed1.fq.gz \
	 --output2  $output_dir/fastq/"$(basename $FQGZ _R1.fastq.gz)".trimmed2.fq.gz \
	 --gzip  --trimns --trimqualities --minlength 20
    echo "Finished trimming of "$FQGZ" "
       
done
