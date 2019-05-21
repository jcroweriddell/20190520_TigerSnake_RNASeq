#!/bin/bash

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=00:05:00
#SBATCH --mem=4GB
#SBATCH -o /home/a1671704/fastdir/tigerSnake_Paper/0_rawData/slurm/gzipfiles_raw_%j.out
#SBATCH -e /home/a1671704/fastdir/tigerSnake_Paper/0_rawData/slurm/gzipfiles_raw_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=a1671704@student.adelaide.edu.au

## set path to directory
RAWDATA=/home/a1671704/fastdir/tigerSnake_Paper/0_rawData/fastq/AGRF_CAGRF17322_H37VKDRXX

## gzip fstq files in directory and run in loop
for ${f1} in ${RAWDATA}/*.fastq
  do
    echo '${f1}'
  done
