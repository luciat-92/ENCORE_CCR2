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
scripts/0_batch_correction_encore.R
```

#### STEP 1:
Data preprocessing of genome-wide single screens (Project score and new screens) and Encore data.
```
scripts/1a_preproc_project_score_raw.R 
scripts/1b_preproc_encore_raw.R
```

#### STEP 2:
Application of [CRISPRcleanR^2](https://github.com/luciat-92/CRISPRcleanRatSquared.git) to each cell line separately. Bash script to run in parallel.
```
sbatch --array=1-32%11 scripts/2_ccr2_run.sh
```

#### STEP 3:
Summarize output across all cell lines and generate relevant plots.
```
scripts/3_summarise.R
```


