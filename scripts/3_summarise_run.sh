#!/bin/bash
#SBATCH --job-name=ccr2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucia.trastulla@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=1-0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/ccr2_summary_%j.out.log
#SBATCH --error=logs/ccr2_summary_%j.err.log
#SBATCH --mem=40G


######################################################################################################################
### Set the environment
######################################################################################################################
module load nlopt/2.7.0-intel-oneapi-mkl-2021.4.0 
module load R/4.1.0

# sbatch scripts/3_summarise_run.sh ORIGINAL c91
# sbatch scripts/3_summarise_run.sh ORIGINAL c92
# sbatch scripts/3_summarise_run.sh ORIGINAL cavg
# sbatch scripts/3_summarise_run.sh BATCH_CORRECTED c91
# sbatch scripts/3_summarise_run.sh BATCH_CORRECTED c92
# sbatch scripts/3_summarise_run.sh BATCH_CORRECTED cavg

GROUP_FOLDER=/group/iorio/lucia/
TYPE_PRE=$1
TYPE_CASNEG=$2

echo "$TYPE_PRE $TYPE_CASNEG" 
mkdir -p ${GROUP_FOLDER}/CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/${TYPE_PRE}/${TYPE_CASNEG}/

Rscript R/3_summarise_run.R \
	--fold_input_dual ${GROUP_FOLDER}/CRISPR_combinatorial/data/encore/DATA_FREEZE_v4_NOV_2023/${TYPE_PRE}/${TYPE_CASNEG}/ \
  --fold_output ${GROUP_FOLDER}/CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/${TYPE_PRE}/${TYPE_CASNEG}/ \
  --fold_CN ${GROUP_FOLDER}datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4_NOV_2023/METADATA_FEB2023/COPY_NUMBER/NEW_COPY_NUMBER/
