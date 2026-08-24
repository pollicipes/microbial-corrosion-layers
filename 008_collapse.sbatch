#!/bin/bash

#SBATCH --output=./log/%x_o.%A_%a
#SBATCH --error=./log/%x_e.%A_%a
#SBATCH -c 4
#SBATCH --mem=6G
#SBATCH --time=06:00:00 
#SBATCH --array=1-11%11

set -Eeuo pipefail;

mkdir -p ./log;
module purge;

# Modify the path accordingly
list="/home/ngp704/data/iron_reducing_bacteria/raw_data/list"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${list}  | awk '{print $1}');

module load fastp;

cd "/home/ngp704/data/iron_reducing_bacteria/raw_data";

fastp -i ${sample}_clean_1.fastq.gz \
      -o ${sample}_clean_1_unmerged.fastq.gz \
      -I ${sample}_clean_2.fastq.gz \
      -O ${sample}_clean_2_unmerged.fastq.gz \
      -m \
      --merged_out ${sample}_clean_collapsed.fastq.gz \
      --overlap_len_require 11
