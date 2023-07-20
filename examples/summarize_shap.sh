#!/usr/bin/env bash
#SBATCH -A naiss2023-5-211
#SBATCH -J nando
#SBATCH -n 5
#SBATCH -t 0-72:00:00
#SBATCH -p core

######################################################
module load bioinfo-tools
module load R_packages/4.2.1 

######################################################
# Fetch one directory from the array based on the task ID (index starts from 0)
echo "Doing file $SLURM_ARRAY_TASK_ID"

Rscript summarize_shap.R

exit 0
######################################################
######################################################
######################################################

# sbatch --array=0-50 summarize_shap.sh		#full speed!
# sbatch --array=0-0 summarize_shap.sh			#testing
