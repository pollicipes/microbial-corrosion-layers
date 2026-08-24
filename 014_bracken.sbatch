#!/bin/bash

#SBATCH --output=./log/%x_o.%A_%a
#SBATCH --error=./log/%x_e.%A_%a
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=2        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=4G
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --array=1-11%11

set -Eeuo pipefail;

module load bracken;

level=$1 # G for Genus; F for Family, etc.

# Modify the path accordingly
list="/home/ngp704/data/iron_reducing_bacteria/list"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${list}  | awk '{print $1}');
DB="/datasets/globe_databases/krakenuniq/standard/20200816/2020-08-16";
INP="/home/ngp704/data/iron_reducing_bacteria/results/krakenuniq_collapsed/"
OUT="/home/ngp704/data/iron_reducing_bacteria/results/bracken"

mkdir -p ${OUT}
cd ${OUT}

# KMER="/datasets/globe_databases/krakenuniq/standard/20200816/2020-08-16/database100mers.kmer_distrib"
# generate_kmer_distribution.py -i ${KMER} -o ${OUT}/database100mers.kmer_distrib

bracken -d ${OUT} -i ${INP}"/"${sample}"_kraken_report_standard" \
	-o ${OUT}"/"${sample}"_bracken_report_"${level}".abundance" \
	-w ${OUT}"/"${sample}"_bracken_report_"${level}".out" \
	-l ${level} # -r  <read length>

# TO COMBINE THE REPORTS:
# (RUN MANUALLY)
combine_bracken_output.py --files ${OUT}/*_bracken_report_${level}.abundance -o bracken_abundances_${level}_ALL
