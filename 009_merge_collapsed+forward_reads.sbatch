#!/bin/bash

#SBATCH --output=./log/%x_o.%A_%a
#SBATCH --error=./log/%x_e.%A_%a
#SBATCH -c 3
#SBATCH --mem=6G
#SBATCH --time=02:00:00 
#SBATCH --array=1-11%11

set -Eeuo pipefail;

mkdir -p ./log;
module purge;

# Modify the path accordingly
list=/home/ngp704/data/iron_reducing_bacteria/list
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${list}  | awk '{print $1}');

cd /home/ngp704/data/iron_reducing_bacteria/raw_data/;

cat ${sample}_clean_collapsed.fastq.gz ${sample}_clean_1_unmerged.fastq.gz > ${sample}_clean_collapsed+forward.fastq.gz;
