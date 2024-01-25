#!/bin/bash
#SBATCH --job-name=ccr2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucia.trastulla@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=1-0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/compare_cas9neg_%j.out.log
#SBATCH --error=logs/compare_cas9neg_%j.err.log
#SBATCH --mem=10G


######################################################################################################################
### Set the environment
######################################################################################################################
module load nlopt/2.7.0-intel-oneapi-mkl-2021.4.0 
module load R/4.1.0

# sbatch scripts/4b_compare_cas9_norm_run.sh ORIGINAL
# sbatch scripts/4b_compare_cas9_norm_run.sh BATCH_CORRECTED

GROUP_FOLDER=/group/iorio/lucia/
TYPE_CORR=$1

echo "$TYPE_CORR" 
COMMON_NAME=${GROUP_FOLDER}/CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/${TYPE_CORR}/
mkdir -p ${COMMON_NAME}/SUMMARY_cas9neg_norm/

Rscript R/4b_compare_cas9_norm_run.R \
  --fold_output ${GROUP_FOLDER}/CRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/${TYPE_CORR}/SUMMARY_cas9neg_norm/ \
  --name_input cavg c91 c92 \
  --fold_input ${COMMON_NAME}/cavg/ALL_CLs/ ${COMMON_NAME}/c91/ALL_CLs/ ${COMMON_NAME}/c92/ALL_CLs/ 