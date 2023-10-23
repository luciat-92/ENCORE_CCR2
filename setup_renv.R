## setup of renv (ccr2 under development)
library("renv")

#### INITIALIZE #####
renv::init(
  settings = list(snapshot.type = "implicit", # default, Capture only packages which appear to be used in your project
                  package.dependency.fields = c("Imports", "Depends", "LinkingTo", "Remotes")), # When explicitly installing a package with install(), what fields should be used to determine that packages dependencies?
)
renv::activate()

#### INSTALL ####
renv::install("devtools")
renv::install("tidyverse")
renv::install("pak")

### IMPORTANT: if not working properly, specify in installation: repos = "https://cran.mirror.garr.it/CRAN/"
# add github credential PAT (only needed if repo is not public)
library(devtools)
library(pak)
gitcreds::gitcreds_set()

pak::pkg_install("luciat-92/CRISPRcleanRatSquared") # if using renv::install, there are problems with the dependencies of ccr

#### SNAPSHOT ####
renv::snapshot()

