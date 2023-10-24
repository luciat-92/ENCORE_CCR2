# Apply systematically across all CLs
library(renv)
renv::activate()
print(.libPaths())

library(CRISPRcleanRatSquared)
library(tidyverse)
library(argparse)
pdf(NULL) # avoid Rplots.pdf from being generated

parser <- ArgumentParser(description = "Run CRISPRcleanR^2 for a specified CL across all available libraries")
parser$add_argument("--CL_name", type = "character", help = "Cell line name as available in CellModelPassport")

args <- parser$parse_args()
CL_name <- args$CL_name

print(paste("CL:", CL_name))

###########################################

source("R/2_auxilary_functions.R")

root_path <- "/group/iorio/lucia/"
input_dual_fold <- sprintf("%sdatasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4/", root_path)
input_single_fold <- sprintf("%sdatasets/PROJECT_SCORE/SINGLE_CL/", root_path)
output_fold <- sprintf("%sCRISPR_combinatorial/", root_path)

# load CL info
model_dual_table <- read_tsv(sprintf("%sBATCH_CORRECTED/model_list.tsv", input_dual_fold), 
                        show_col_types = FALSE) %>%
  dplyr::filter(model_name_CMP == CL_name)

model_single_table <- read_tsv(sprintf("%smodel_list.tsv", input_single_fold), 
                               show_col_types = FALSE) %>%
  dplyr::filter(model_id_CMP %in% model_dual_table$model_id_CMP) %>% # only CLs in dual
  dplyr::distinct(model_name_CMP_library, .keep_all = TRUE) %>% 
  dplyr::select(model_name_CMP_library, model_name_CMP, library, model_id_CMP)

if (nrow(model_single_table) == 0) { 
  stop("No single screen available")
}

model_dual_table <- model_dual_table %>%
  dplyr::mutate(model_name_CMP_lib = paste(model_name_CMP, lib, sep = "_")) %>%
  dplyr::distinct(model_name_CMP_lib, .keep_all = TRUE) %>%
  dplyr::select(model_name_CMP_lib, model_name_CMP, lib, model_id_CMP)

##############
#### COLO ####
##############

model_dual_COLO <- model_dual_table[grepl("COLO", model_dual_table$lib),]
if (nrow(model_dual_COLO) > 0) {
  
  model_single_COLO <- model_single_table %>% 
    dplyr::filter(model_id_CMP %in% model_dual_COLO$model_id_CMP)
  
  COLO_input_list <- create_input_list(
    model_dual = model_dual_COLO,
    model_single = model_single_COLO,
    tissue = "COLO", 
    copy_number_file = "METADATA_FEB2023/COPY_NUMBER/NEW_COPY_NUMBER/MERGED_SEGMENT_COPYNUMBER.txt",
    input_fold = input_dual_fold)
   
  # get single
  fn_single <- COLO_input_list$single$count_file
  library_single <- get(base::load(COLO_input_list$single$library_file))
  # get dual
  input_dual <- get_input_data.v1(param_list = COLO_input_list$dual)
  CNA <- input_dual$CNA
  dual_logFC <- input_dual$result
  dual_library <- input_dual$library
  CL_name <- input_dual$CL_name
  out_fold <- input_dual$out_fold
  
  if (any(CNA$CN_category == "Amp")) {
    CN_thr <- round(min(CNA$C[CNA$CN_category == "Amp"]))  
  }else{
    CN_thr <- 8
  }
    
  print(paste0("###### ", CL_name, " (COLO) ", " ######"))
  
  system(sprintf("mkdir -p %s/COLO/",  out_fold))
    
  output <- NULL
  try(output <- ccr2.run_complete(
    filename_single = fn_single,
    EXPname = CL_name,
    libraryAnnotation_single = library_single,
    display = FALSE,
    outdir = sprintf("%s/COLO/", out_fold),
    saveToFig = TRUE,
    libraryAnnotation_dual = dual_library,
    dual_logFC = dual_logFC,
    correctGW = NULL,
    CNA = CNA,
    CN_thr = CN_thr,
    saveFormat = "png",
    split_graph = TRUE)
  )
  # save output
  save(output, file = sprintf("%s/COLO/%s_dual_correctedFCs.RData",  out_fold, CL_name))
  
}

##############
#### BRCA ####
##############

model_dual_BRCA <- model_dual_table[grepl("BRCA", model_dual_table$lib),]
if (nrow(model_dual_BRCA) > 0) {
  
  model_single_BRCA <- model_single_table %>% 
    dplyr::filter(model_id_CMP %in% model_dual_BRCA$model_id_CMP)
  
  BRCA_input_list <- create_input_list(
    model_dual = model_dual_BRCA,
    model_single = model_single_BRCA,
    tissue = "BRCA", 
    copy_number_file = "METADATA_FEB2023/COPY_NUMBER/NEW_COPY_NUMBER/MERGED_SEGMENT_COPYNUMBER.txt",
    input_fold = input_dual_fold)
  
  # get single
  fn_single <- BRCA_input_list$single$count_file
  library_single <- get(base::load(BRCA_input_list$single$library_file))
  # get dual
  input_dual <- get_input_data.v1(param_list = BRCA_input_list$dual)
  CNA <- input_dual$CNA
  dual_logFC <- input_dual$result
  dual_library <- input_dual$library
  CL_name <- input_dual$CL_name
  out_fold <- input_dual$out_fold
  
  if (any(CNA$CN_category == "Amp")) {
    CN_thr <- round(min(CNA$C[CNA$CN_category == "Amp"]))  
  }else{
    CN_thr <- 8
  }
  
  print(paste0("###### ", CL_name, " (BRCA) ", " ######"))
  
  system(sprintf("mkdir -p %s/BRCA/",  out_fold))
  
  output <- NULL
  try(output <- ccr2.run_complete(
    filename_single = fn_single,
    EXPname = CL_name,
    libraryAnnotation_single = library_single,
    display = FALSE,
    outdir = sprintf("%s/BRCA/", out_fold),
    saveToFig = TRUE,
    libraryAnnotation_dual = dual_library,
    dual_logFC = dual_logFC,
    correctGW = NULL,
    CNA = CNA,
    CN_thr = CN_thr,
    saveFormat = "png",
    split_graph = TRUE)
  )
  # save output
  save(output, file = sprintf("%s/BRCA/%s_dual_correctedFCs.RData",  out_fold, CL_name))
  
}






