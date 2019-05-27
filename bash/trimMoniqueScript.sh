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
module load AdapterRemoval

## Script for FastQC and adapter removal
working=$(pwd)

data_dir=/fast/users/a1662801/20190520_TigerSnake_RNASeq/0_rawData/fastq/
output_dir=/20190520_TigerSnake_RNASeq/1_trimmedData/fastq

# Run FastQC on each read pair
#echo "Starting fastqc"

#fastqc -o $data_dir/FastQC $data_dir/fastq/*.fastq.gz 

#echo "Finished fastqc"

for f in ${data_dir}/*R1.fastq
  do
    f2=${f%R1.fastq}R2.fastq
    echo "found file ${f}"
    echo "how about ${f2}"
    AdapterRemoval \
    --file1 ${f} \
    --file2 ${f2} \
    --output1 ${output_dir}/$(basename ${f}) \
    --output2 ${output_dir}/$(basename ${f2}) \
    --threads 16 \
    --trimqualities \
    --minquality 30 \
    --trimns \
    --minlength 35 
done
