# script to convert data in a separate file for each cell line
library(tidyverse)
setwd("/group/iorio/lucia/")

############
# data have been batch corrected and all libraries are combined
fold_data <- "CRISPR_combinatorial/data/encore/DATA_FREEZE_v4_NOV_2023/"
file_CMP_ann <- "datasets/CMP_PORTAL/model_annotation/model_list_20230923.csv"
fold_out <- "CRISPR_combinatorial/data/encore/DATA_FREEZE_v4_NOV_2023/"
fold_CN <- "datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4_NOV_2023/METADATA_FEB2023/"
fold_input_single <- "datasets/PROJECT_SCORE/SINGLE_CL/"
libs_name <- c("COLO", "BRCA")
cneg_label <- c("c91", "c92", "cavg")
###########

# get cell model passport
CMP_table <- read_csv(file_CMP_ann) %>%
  dplyr::select(model_id, sample_id, model_name, synonyms, tissue, cancer_type, 
                tissue_status, COSMIC_ID, BROAD_ID, CCLE_ID) %>%
  dplyr::rename(model_id_CMP = model_id, 
                sample_id_CMP = sample_id)

# combine copy number in a unique file
COLO_CNA <- read_table(sprintf("%sCOPY_NUMBER/NEW_COPY_NUMBER/COLO/COLO_MERGED_SEGMENT_COPYNUMBER.txt", fold_CN))
BRCA_CNA <- read_table(sprintf("%sCOPY_NUMBER/NEW_COPY_NUMBER/BRCA/BRCA_MERGED_SEGMENT_COPYNUMBER.txt", fold_CN))

tot_CNA <- bind_rows(COLO_CNA, BRCA_CNA)
write_tsv(x = tot_CNA, 
          file = sprintf("%sCOPY_NUMBER/NEW_COPY_NUMBER/MERGED_SEGMENT_COPYNUMBER.txt", fold_CN), 
          col_names = TRUE)

# load single GW CLs
model_single_table <- read_tsv(sprintf("%smodel_list.tsv", fold_input_single), 
                               show_col_types = FALSE) 

# get sample list encore
### batch corrected ####
for (idx_cneg in seq_len(length(cneg_label))) {
  
  print(paste0("##### ", cneg_label[idx_cneg], " #####"))
  input <- list()
  samples_name <- list()
  
  for (idx in seq_len(length(libs_name))) {
    
    lib_id <- libs_name[idx]
    print(paste0("##### ", lib_id, " #####"))
    file_name <- sprintf("%s/BATCH_CORRECTED/%s/%s_FINAL_EXACT_logFC_sgRNA.txt", 
                         fold_data, cneg_label[idx_cneg], lib_id)
    
    input[[idx]] <- readr::read_tsv(file_name, 
                                    col_types = readr::cols(.default = "?",
                                                            sgRNA1_WGE_ID = "c", 
                                                            sgRNA2_WGE_ID = "c"), 
                                    show_col_types = FALSE)
    
    samples_name[[idx]] <- data.frame(id = colnames(input[[idx]])) %>%
      dplyr::filter(!grepl("sgRNA", id) & !id %in% c("ID", "lib", 
                                                     "ID_lib", "Note1", 
                                                     "Note2", "Gene_Pair", 
                                                     "Gene1", "Gene2", 
                                                     "SEQ_pair")) %>%
      dplyr::mutate(lib = lib_id, model_name = id)
    
  }
  samples_name <- bind_rows(samples_name)
  
  model_dual_table <- left_join(samples_name, CMP_table, by = "model_name") %>%
    mutate(model_name_uppercase = str_replace_all(model_name, "[-]", "")) %>%
    mutate(model_name_uppercase = toupper(model_name_uppercase)) %>% 
    dplyr::rename(model_name_CMP = model_name)
  
  # save
  write_tsv(x = model_dual_table,
            file = sprintf("%s/BATCH_CORRECTED/%s/model_list.tsv", fold_out, cneg_label[idx_cneg]), 
            col_names = T)
  
  ## get list of CLs having matching single screen info ##
  model_dual_table <- model_dual_table %>%
    dplyr::mutate(model_name_CMP_lib = paste(model_name_CMP, lib, sep = "_")) %>%
    dplyr::distinct(model_name_CMP_lib, .keep_all = TRUE) %>%
    dplyr::select(model_name_CMP_lib, model_name_CMP, lib, model_id_CMP)
  
  model_single_table_tmp <- model_single_table %>%
    dplyr::filter(model_id_CMP %in% model_dual_table$model_id_CMP) %>% # only CLs in encore
    dplyr::distinct(model_name_CMP_library, .keep_all = TRUE) %>% 
    dplyr::select(model_name_CMP_library, model_name_CMP, library, model_id_CMP)
  
  # create summary CL
  CL_match_summary <- model_dual_table %>%
    dplyr::distinct(model_name_CMP, .keep_all = TRUE) %>%
    dplyr::mutate(single_screen = model_name_CMP %in% model_single_table_tmp$model_name_CMP) %>%
    dplyr::select(-model_name_CMP_lib, -lib)
  
  write_tsv(sprintf("%s/BATCH_CORRECTED/%s/models_match_single_screen.tsv", fold_out, cneg_label[idx_cneg]), 
            x = CL_match_summary, 
            col_names = TRUE)
  
}


### original ###
for (idx_cneg in seq_len(length(cneg_label))) {
  
  print(paste0("##### ", cneg_label[idx_cneg], " #####"))
  input <- list()
  samples_name <- list()
  
  for (idx in seq_len(length(libs_name))) {
    
    lib_id <- libs_name[idx]
    print(paste0("##### ", lib_id, " #####"))
    file_name <- sprintf("%s/ORIGINAL/%s/%s_FINAL_EXACT_logFC_sgRNA.txt", 
                         fold_data, cneg_label[idx_cneg], lib_id)
    
    input[[idx]] <- readr::read_tsv(file_name, 
                                    col_types = readr::cols(.default = "?",
                                                            sgRNA1_WGE_ID = "c", 
                                                            sgRNA2_WGE_ID = "c"), 
                                    show_col_types = FALSE)
    
    samples_name[[idx]] <- data.frame(id = colnames(input[[idx]])) %>%
      dplyr::filter(!grepl("sgRNA", id) & !id %in% c("ID", "lib", 
                                                     "ID_lib", "Note1", 
                                                     "Note2", "Gene_Pair", 
                                                     "Gene1", "Gene2", 
                                                     "SEQ_pair")) %>%
      dplyr::mutate(lib = lib_id, model_name = id)
    
  }
  samples_name <- bind_rows(samples_name)
  
  model_dual_table <- left_join(samples_name, CMP_table, by = "model_name") %>%
    mutate(model_name_uppercase = str_replace_all(model_name, "[-]", "")) %>%
    mutate(model_name_uppercase = toupper(model_name_uppercase)) %>% 
    dplyr::rename(model_name_CMP = model_name)
  
  # save
  write_tsv(x = model_dual_table,
            file = sprintf("%s/ORIGINAL/%s/model_list.tsv", fold_out, cneg_label[idx_cneg]), 
            col_names = T)
  
  ## get list of CLs having matching single screen info ##
  model_dual_table <- model_dual_table %>%
    dplyr::mutate(model_name_CMP_lib = paste(model_name_CMP, lib, sep = "_")) %>%
    dplyr::distinct(model_name_CMP_lib, .keep_all = TRUE) %>%
    dplyr::select(model_name_CMP_lib, model_name_CMP, lib, model_id_CMP)
  
  model_single_table_tmp <- model_single_table %>%
    dplyr::filter(model_id_CMP %in% model_dual_table$model_id_CMP) %>% # only CLs in encore
    dplyr::distinct(model_name_CMP_library, .keep_all = TRUE) %>% 
    dplyr::select(model_name_CMP_library, model_name_CMP, library, model_id_CMP)
  
  # create summary CL
  CL_match_summary <- model_dual_table %>%
    dplyr::distinct(model_name_CMP, .keep_all = TRUE) %>%
    dplyr::mutate(single_screen = model_name_CMP %in% model_single_table_tmp$model_name_CMP) %>%
    dplyr::select(-model_name_CMP_lib, -lib)
  
  write_tsv(sprintf("%s/ORIGINAL/%s/models_match_single_screen.tsv", fold_out, cneg_label[idx_cneg]), 
            x = CL_match_summary, 
            col_names = TRUE)
  
}







