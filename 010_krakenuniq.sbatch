#!/bin/bash

#SBATCH --output=./log/%x_o.%A_%a
#SBATCH --error=./log/%x_e.%A_%a
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=20        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=300G
#SBATCH --time=120:00:00          # total run time limit (HH:MM:SS)
#SBATCH --array=1,2,5,6,7,8,9,10,11%9

module load krakenuniq/1.0.4;
set -Eeuo pipefail;

# Modify the path accordingly
list="/home/ngp704/data/iron_reducing_bacteria/list"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${list}  | awk '{print $1}');

data="/home/ngp704/data/iron_reducing_bacteria/raw_data/";
cd ${data};

DB="/datasets/globe_databases/krakenuniq/standard/20200816/2020-08-16";
tipo="standard";

OUT="/home/ngp704/data/iron_reducing_bacteria/results/krakenuniq_collapsed/";
mkdir -p ${OUT};

krakenuniq --report-file ${OUT}"/"${sample}"_kraken_report_"${tipo} --threads 20 --db ${DB} ${sample}"_clean_collapsed+forward.fastq.gz";

