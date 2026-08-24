#!/bin/bash

#SBATCH --output=./log/%x_o.%A_%a
#SBATCH --error=./log/%x_e.%A_%a
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --time=02:00:00 
#SBATCH --array=1-11%11

set -Eeuo pipefail;

mkdir -p ./log;
module purge;
module load metaphlan;
module load bowtie2;
module load samtools;

# Modify the path accordingly
list="/home/ngp704/data/iron_reducing_bacteria/list"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${list}  | awk '{print $1}');

WD="/home/ngp704/data/iron_reducing_bacteria/metaphlan/out/"
mkdir -p ${WD};

FQ="/home/ngp704/data/iron_reducing_bacteria/raw_data/"

# metaphlan ${FQ}/${sample}_clean_collapsed+forward.fastq.gz \
# 	  --bowtie2db /datasets/globe_databases/metaphlan/20241118/ \
# 	  --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
# 	  --bowtie2out ${WD}/${sample}.outfmt.bz2 \
# 	  --nproc 8 \
# 	  --input_type fastq \
# 	  -o ${WD}/${sample}.profile.txt


metaphlan ${WD}/${sample}.outfmt.bz2 \
 	  --bowtie2db /datasets/globe_databases/metaphlan/20241118/ \
 	  --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
	  --input_type bowtie2out \
	  -t rel_ab_w_read_stats \
	  --nproc 4 > ${WD}/${sample}.reads_profile.txt;

# samtools view -bS ${WD}/${sample}.sam > ${WD}/${sample}.bam; 
