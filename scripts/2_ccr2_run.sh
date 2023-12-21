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

# sbatch --array=1-33%11 scripts/2_ccr2_run.sh ORIGINAL c91
# sbatch --array=1-33%11 scripts/2_ccr2_run.sh ORIGINAL c92
# sbatch --array=1-33%11 scripts/2_ccr2_run.sh ORIGINAL cavg
# sbatch --array=1-33%11 scripts/2_ccr2_run.sh BATCH_CORRECTED c91
# sbatch --array=1-33%11 scripts/2_ccr2_run.sh BATCH_CORRECTED c92
# sbatch --array=1-33%11 scripts/2_ccr2_run.sh BATCH_CORRECTED cavg

GROUP_FOLDER=/group/iorio/lucia/
idx=${SLURM_ARRAY_TASK_ID}
TYPE_PRE=$1
TYPE_CASNEG=$2

CLs_info="${GROUP_FOLDER}/CRISPR_combinatorial/data/encore/DATA_FREEZE_v4_NOV_2023/${TYPE_PRE}/${TYPE_CASNEG}/models_match_single_screen.tsv"
CLs=($(tail -n +2 "$CLs_info" | awk -F"\t" '{print $1}'))
CL_name=$(eval echo "\${CLs[${idx}-1]}")

echo "$TYPE_PRE $TYPE_CASNEG $CL_name" 
mkdir -p ${GROUP_FOLDER}/CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/${TYPE_PRE}/${TYPE_CASNEG}/
  
Rscript scripts/2_ccr2_perCL_run.R \
	--CL_name ${CL_name} \
	--fold_input_dual ${GROUP_FOLDER}/CRISPR_combinatorial/data/encore/DATA_FREEZE_v4_NOV_2023/${TYPE_PRE}/${TYPE_CASNEG}/ \
  --fold_input_single ${GROUP_FOLDER}/datasets/PROJECT_SCORE/SINGLE_CL/ \
  --fold_output ${GROUP_FOLDER}/CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/${TYPE_PRE}/${TYPE_CASNEG}/ \
  --fold_CN ${GROUP_FOLDER}datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4_NOV_2023/METADATA_FEB2023/COPY_NUMBER/NEW_COPY_NUMBER/
