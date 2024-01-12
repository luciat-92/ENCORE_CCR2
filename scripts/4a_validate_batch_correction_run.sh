#!/bin/bash
#SBATCH --job-name=ccr2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucia.trastulla@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=1-0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/validate_batchcorr_%j.out.log
#SBATCH --error=logs/validate_batchcorr_%j.err.log
#SBATCH --mem=50G


######################################################################################################################
### Set the environment
######################################################################################################################
module load nlopt/2.7.0-intel-oneapi-mkl-2021.4.0 
module load R/4.1.0

# sbatch scripts/4a_validate_batch_correction_run.sh c91
# sbatch scripts/4a_validate_batch_correction_run.sh c92
# sbatch scripts/4a_validate_batch_correction_run.sh cavg

GROUP_FOLDER=/group/iorio/lucia/
TYPE_CASNEG=$1

echo "$TYPE_CASNEG" 
mkdir -p ${GROUP_FOLDER}/CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/SUMMARY_ORIGINAL_BATCH_CORRECTED/${TYPE_CASNEG}/

Rscript R/4a_validate_batch_correction_run.R \
  --fold_output ${GROUP_FOLDER}/CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/SUMMARY_ORIGINAL_BATCH_CORRECTED/${TYPE_CASNEG}/ \
  --fold_original ${GROUP_FOLDER}CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/ORIGINAL/${TYPE_CASNEG}/ALL_CLs/ \
  --fold_batchcorr ${GROUP_FOLDER}CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/BATCH_CORRECTED/${TYPE_CASNEG}/ALL_CLs/
