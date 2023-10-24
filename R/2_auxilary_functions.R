### auxiliary functions ###
create_input_list <- function(model_dual,
                              model_single,
                              tissue, 
                              input_fold, 
                              copy_number_file){
  
  CL_name <- model_dual$model_name_CMP
  # get dual
  if (length(CL_name) > 1) {stop("Only one value per cell line is permitted")}
  
  result_file <- sprintf("BATCH_CORRECTED/%s_FINAL_EXACT_logFC_sgRNA_ComBatCorrectionLIBs.txt", tissue)
  names(result_file) <- "result_file"
  if (tissue == "COLO") {
    library_file <- sprintf("BATCH_CORRECTED/ENCORE_GI_COREAD_Library_ALL.txt")
  }
  if (tissue == "BRCA") {
    library_file <- sprintf("BATCH_CORRECTED/ENCORE_GI_BRCA_Library_ALL.txt")
  }
  names(library_file) <- "library_file"
  
  out_fold <- sprintf("%sCRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4/%s/", root_path, CL_name)
  names(out_fold) <- "out_fold"
  names(input_fold) <- "input_fold"
  names(copy_number_file) <- "copy_number_file"
  names(CL_name) <- "CL_name"
  
  dual_input <- as.list(c(input_fold, 
                          copy_number_file, 
                          result_file, 
                          library_file, 
                          CL_name, 
                          out_fold))
  
  # get single
  CL_name <- model_dual$model_name_CMP
  if (CL_name %in% model_single$model_name_CMP) {
    tmp <- model_single %>% filter(model_name_CMP %in% CL_name)
    if (nrow(tmp) > 1) {
      tmp <- tmp %>% filter(library == "V1.1")
    }
    
    library_file <- sprintf("%sCRISPR_libraries/data/KY_Library_%s_hg38.RData", root_path, tolower(tmp$library))
    names(library_file) <- "library_file"
    count_file <- sprintf("%sdatasets/PROJECT_SCORE/SINGLE_CL/%s_counts.tsv", root_path, tmp$model_name_CMP_library)
    names(count_file) <- "count_file"
    names(CL_name) <- "CL_name"
    single_input <- as.list(c(count_file, 
                              library_file, 
                              CL_name))
    
  }else{
    single_input <- NULL
  }
  return(list(dual = dual_input, single = single_input))
}

