#!/bin/bash

#SBATCH --output=./log/%x_o.%A_%a
#SBATCH --error=./log/%x_e.%A_%a
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=02:00:00 
#SBATCH --array=1-1%1
### SBATCH --begin=2023-07-22T04:30:00

set -Eeuo pipefail;

mkdir -p ./log;
module purge;
module load fastp;
origin=$1

# Modify the path accordingly
list="/projects/mjolnir1/people/fwc915/Iron_reducing_"${origin}"/raw_data/sample_list.txt";
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${list}  | awk '{print $1}');

RAW="/projects/mjolnir1/people/fwc915/Iron_reducing_"${origin}"/raw_data/";
CLEAN="/home/ngp704/data/iron_reducing_bacteria/raw_data";
FASTQC="/home/ngp704/data/iron_reducing_bacteria/raw_data/fastqc"

mkdir -p ${CLEAN};
mkdir -p ${FASTQC};

fastp -i ${RAW}"/"${sample}"_1.fq.gz" \
      -o ${CLEAN}"/"${sample}"_clean_1.fq.gz" \
      -I ${RAW}"/"${sample}"_2.fq.gz" \
      -O ${CLEAN}"/"${sample}"_clean_2.fq.gz" \
      --detect_adapter_for_pe \
      -w 8 \
      -f 10 \
      -l 30 \
      -h ${CLEAN}"/"${sample}".html"

# # FASTQC ####
module load openjdk perl fastqc
cd ${FASTQC};
fastqc -t 8 ${RAW}"/"${sample}"_1.fq.gz" -o ${FASTQC};
fastqc -t 8 ${RAW}"/"${sample}"_2.fq.gz" -o ${FASTQC};

fastqc -t 8 ${CLEAN}"/"${sample}"_clean_1.fq.gz" -o ${FASTQC};
fastqc -t 8 ${CLEAN}"/"${sample}"_clean_2.fq.gz" -o ${FASTQC};
