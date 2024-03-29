# ENCORE_CCR2
Application of CRISPRcleanR^2 to ENCORE dual combinatorial screens.
- Data Release: DATA_FREEZE_v4_Nov2023

## SPECIFICS:
To reproduce results use R/4.1.0 and renv 0.15.5 package
To install libraries:
```R
library(renv)
renv::restore()
```

## WORKFLOW

#### STEP 0: 
Batch correction of ENCORE data via ComBat. 
**ATM: outputs both not corrected and corrected format. Repeated per c91, c92 and cavg versions**.
The output is a table with one guide pair in each library per row and CLs + guide pair info per columns. The data is in the logFC format, replicates for each CL are collapsed taking the mean.
```
sbatch scripts/0_batch_correction_run.sh
```

#### STEP 1:
Data preprocessing of genome-wide single screens (Project score and new screens) and Encore data.
```
sbatch scripts/1_preproc_run.sh 
```

#### STEP 2:
Application of [CRISPRcleanR^2](https://github.com/luciat-92/CRISPRcleanRatSquared.git) to each cell line, type of preprocessing (ORIGINAL or BATCH CORRECTED) and type of cas9Negative, separately. Bash script to run in parallel.
```
sbatch --array=1-32 scripts/2_ccr2_run.sh ORIGINAL c91
sbatch --array=1-32 scripts/2_ccr2_run.sh ORIGINAL c92
sbatch --array=1-32 scripts/2_ccr2_run.sh ORIGINAL cavg
sbatch --array=1-32 scripts/2_ccr2_run.sh BATCH_CORRECTED c91
sbatch --array=1-32 scripts/2_ccr2_run.sh BATCH_CORRECTED c92
sbatch --array=1-32 scripts/2_ccr2_run.sh BATCH_CORRECTED cavg
```

#### STEP 3:
Summarize output across all cell lines and generate relevant plots.
```
sbatch scripts/3_summarise_run.sh ORIGINAL c91
sbatch scripts/3_summarise_run.sh ORIGINAL c92
sbatch scripts/3_summarise_run.sh ORIGINAL cavg
sbatch scripts/3_summarise_run.sh BATCH_CORRECTED c91
sbatch scripts/3_summarise_run.sh BATCH_CORRECTED c92
sbatch scripts/3_summarise_run.sh BATCH_CORRECTED cavg
```

#### STEP 4:
Investigate results: ComBat batch correction VS Original; comparison cavg, c91, c92
```
sbatch scripts/4a_validate_batch_correction_run.sh c91
sbatch scripts/4a_validate_batch_correction_run.sh c92
sbatch scripts/4a_validate_batch_correction_run.sh cavg

sbatch scripts/4b_compare_cas9_norm_run.sh ORIGINAL
sbatch scripts/4b_compare_cas9_norm_run.sh BATCH_CORRECTED
```