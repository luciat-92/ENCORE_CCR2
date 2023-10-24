#!/bin/bash
#SBATCH --job-name=ccr2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucia.trastulla@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=1-0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/ccr2_%A_CL_idx%a.out.log
#SBATCH --error=logs/ccr2_%A_CL_idx%a.err.log
#SBATCH --mem=10G


######################################################################################################################
### Set the environment
######################################################################################################################
module load nlopt/2.7.0-intel-oneapi-mkl-2021.4.0 
module load R/4.1.0

# sbatch --array=1-33%11 scripts/2_ccr2_run.sh

GROUP_FOLDER=/group/iorio/lucia/

idx=${SLURM_ARRAY_TASK_ID}
CLs_info="${GROUP_FOLDER}/datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4/BATCH_CORRECTED/models_match_single_screen.tsv"
CLs=($(tail -n +2 "$CLs_info" | awk -F"\t" '{print $1}'))
CL_name=$(eval echo "\${CLs[${idx}-1]}")

Rscript scripts/2_ccr2_perCL_run.R \
	--CL_name ${CL_name} \
