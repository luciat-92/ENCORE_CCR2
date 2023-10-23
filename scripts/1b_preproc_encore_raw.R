# script to convert data in a separate file for each cell line
library(tidyverse)

setwd("/group/iorio/lucia/")

# data have been batch corrected and all libraries are combined
fold_data <- "datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4/"
file_CMP_ann <- "datasets/CMP_PORTAL/model_annotation/model_list_20230923.csv"
fold_out <- "datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4/BATCH_CORRECTED/ANALYSIS/"
libs_name <- c("COLO", "BRCA")

input <- list()
samples_name <- list()
for (idx in seq_len(length(libs_name))) {
  
  print(idx)
  lib_id <- libs_name[idx]
  file_name <- sprintf("%s/BATCH_CORRECTED/%s_FINAL_EXACT_logFC_sgRNA_ComBatCorrectionLIBs.txt", 
                       fold_data, lib_id)
  
  input[[idx]] <- readr::read_tsv(file_name, 
                                 col_types = readr::cols(.default = "?",
                                                         sgRNA1_WGE_ID = "c", 
                                                         sgRNA2_WGE_ID = "c"), 
                                 show_col_types = FALSE)
  
  samples_name[[idx]] <- data.frame(id = colnames(input[[idx]])) %>%
    dplyr::filter(!grepl("sgRNA", id) & !id %in% c("ID", "lib", 
                                                   "ID_lib", "Note1", 
                                                   "Note2", "Gene_pair", 
                                                   "Gene1", "Gene2", 
                                                   "SEQ_pair")) %>%
    dplyr::mutate(lib = lib_id, model_name = id)
  
}
samples_name <- bind_rows(samples_name)

# get cell model passport
CMP_table <- read_csv(file_CMP_ann) %>%
  dplyr::select(model_id, sample_id, model_name, synonyms, tissue, cancer_type, 
                tissue_status, COSMIC_ID, BROAD_ID, CCLE_ID) %>%
  dplyr::rename(model_id_CMP = model_id, 
                sample_id_CMP = sample_id)

model_dual_table <- left_join(samples_name, CMP_table, by = "model_name") %>%
  mutate(model_name_uppercase = str_replace_all(model_name, "[-]", "")) %>%
  mutate(model_name_uppercase = toupper(model_name_uppercase))

# save
write_tsv(x = model_dual_table,
          file = sprintf("%smodel_list.tsv", fold_out), 
          col_names = T)

# combine copy number in a unique file
COLO_CNA <- read_table(sprintf("%sMETADATA_FEB2023/COPY_NUMBER/NEW_COPY_NUMBER/COLO/COLO_MERGED_SEGMENT_COPYNUMBER.txt", fold_data))
BRCA_CNA <- read_table(sprintf("%sMETADATA_FEB2023/COPY_NUMBER/NEW_COPY_NUMBER/BRCA/BRCA_MERGED_SEGMENT_COPYNUMBER.txt", fold_data))

tot_CNA <- bind_rows(COLO_CNA, BRCA_CNA)
write_tsv(x = tot_CNA, 
          file = sprintf("%sMETADATA_FEB2023/COPY_NUMBER/NEW_COPY_NUMBER/MERGED_SEGMENT_COPYNUMBER.txt", fold_data), 
                         col_names = TRUE)

## get list of CLs having matching single screen info ##
input_single_fold <- "datasets/PROJECT_SCORE/SINGLE_CL/"
model_single_table <- read_tsv(sprintf("%smodel_list.tsv", input_single_fold), 
                              show_col_types = FALSE) %>%
  dplyr::filter(model_id_CMP %in% model_dual_table$model_id_CMP) %>% # only CLs in encore
  dplyr::distinct(model_name_CMP_library, .keep_all = TRUE) %>% 
  dplyr::select(model_name_CMP_library, model_name_CMP, library, model_id_CMP)

model_dual_table <- model_dual_table %>%
  dplyr::mutate(model_name_CMP_lib = paste(model_name, lib, sep = "_")) %>%
  dplyr::distinct(model_name_CMP_lib, .keep_all = TRUE) %>%
  dplyr::select(model_name_CMP_lib, model_name, lib, model_id_CMP)

# create summary CL
CL_match_summary <- model_dual_table %>%
  dplyr::distinct(model_name, .keep_all = TRUE) %>%
  dplyr::mutate(single_screen = model_name %in% model_single_table$model_name_CMP) %>%
  dplyr::select(-model_name_CMP_lib, -lib)

write_tsv("datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4/BATCH_CORRECTED/models_match_single_screen.tsv", 
          x = CL_match_summary, col_names = TRUE)

