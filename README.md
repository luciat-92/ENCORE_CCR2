# ENCORE_CCR2
Application of CRISPRcleanR^2 to ENCORE dual combinatorial screens.
- Data Release: DATA_FREEZE_v4

## SPECIFICS:
To reproduce results use R/4.1.0 and renv 0.15.5 package
To install libraries:
```R
library(renv)
renv::restore()
```

## WORKFLOW

#### STEP 1:
Data preprocessing of genome-wide single screens (Project score and new screens) and Encore data.
Encore data was previously batch corrected via [ENCORE_batch_correction github](https://github.com/luciat-92/ENCORE_batch_correction.git). 
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


