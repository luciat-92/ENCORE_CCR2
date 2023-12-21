#!/bin/bash
#SBATCH --job-name=batch_corr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucia.trastulla@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/batch_correction_%j.out.log
#SBATCH --error=logs/batch_correction_%j.err.log
#SBATCH --mem=40G

######################################################################################################################
### Set the environment
######################################################################################################################
module load nlopt/2.7.0-intel-oneapi-mkl-2021.4.0 
module load R/4.1.0

Rscript scripts/0_batch_correction_encore.R