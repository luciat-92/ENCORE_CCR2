## setup of renv (ccr2 under development)
library("renv")

#### INITIALIZE #####
#### REDO AND SPECIFY THE CRAN REPO!!!! ####

renv::init(
  settings = list(snapshot.type = "implicit", # default, Capture only packages which appear to be used in your project
                  package.dependency.fields = c("Imports", "Depends", "LinkingTo", "Remotes")), # When explicitly installing a package with install(), what fields should be used to determine that packages dependencies?
)
renv::activate()

#### INSTALL ####
options(repos = BiocManager::repositories())
renv::install("devtools")
renv::install("tidyverse")
renv::install("pak")
renv::install("argparse")
renv::install("pheatmap")
renv::install("R.utils")
renv::install("effectsize")
# biocoductor
renv::install("bioc::Biobase")
renv::install("sva")

# add github credential PAT (only needed if repo is not public)
library(devtools)
library(pak)
gitcreds::gitcreds_set()

pak::pkg_install("luciat-92/CRISPRcleanRatSquared") # if using renv::install, there are problems with the dependencies of ccr

#### SNAPSHOT ####
renv::snapshot()

