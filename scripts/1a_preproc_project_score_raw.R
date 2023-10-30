# convert project score output to the proper format
library(tidyverse)
library(readxl)
library(stringr)
library(dplyr)

setwd("/group/iorio/lucia/")

fold <- "datasets/PROJECT_SCORE/"
# get supplementary table
CL_pub_summary <- readxl::read_xlsx(sprintf("%sRELEASE_1/Behan_Nature_2019/Supplementary Table 1 - CellLines_dataAv_screeDet_QCscores.xlsx", fold), skip = 2, trim_ws = TRUE, n_max = 339)
CL_pub_summary <- CL_pub_summary[, c("Cell\r\nModel Passports\r\nIdentifier", 
                                      "Name",
                                      "Experiment Identifier", 
                                      "sgRNA library", 
                                      "N of Replicates",
                                      "Replicate\r\nPassing\r\nLow Lev QC", 
                                      "Passing High Level Quality Control Assessment \r\n(Analysis Set)")]

df <- list()
for (id_row in seq_len(nrow(CL_pub_summary))) {
  tmp <- CL_pub_summary[id_row, , drop = FALSE]
  df[[id_row]] <- data.frame(
    model_name = rep(pull(tmp[,2]), pull(tmp[,5])), 
    model_id_CMP = rep(pull(tmp[,1]), pull(tmp[,5])), 
    complete_name = str_trim(str_split_fixed(tmp[,3], pattern = "[,]", n = Inf)), 
    passed_low_QC = as.logical(str_trim(str_split_fixed(tmp[,6], pattern = "[,]", n = Inf))), 
    passed_high_QC = rep(pull(tmp[,7]), pull(tmp[,5])),
    library = str_trim(str_split_fixed(tmp[,4], pattern = "[,]", n = Inf)))
}
CL_pub_summary <- bind_rows(df) %>%
  filter(passed_high_QC & (passed_low_QC | is.na(passed_low_QC)))

# get raw counts
# NOTE: there are complete replicates (H716), not in the list, why? 
files <- dir(fold, 
             recursive = TRUE, 
             full.names = FALSE, 
             pattern = "\\.read_count.tsv.gz$")

CL_data_meta <- data.frame(
  location = files, 
  complete_name = str_extract(string = files, pattern = "(?<=/)[^/]+(?=\\.read_count\\.tsv\\.gz)")
  ) %>%
  mutate(model_name_file = str_split_i(string = complete_name, pattern = "_", i = 1))
CL_data_meta$rep <-  sapply(str_split(string = CL_data_meta$complete_name, pattern = "_"), 
         function(x) paste0(x[-1], collapse = "_"))

# filter for those in pub summary
CL_data_meta_QC <- inner_join(CL_data_meta, 
                              CL_pub_summary, 
                              by = "complete_name") %>%
  mutate(origin_data = "PROJECT_SCORE_RELEASE_1")

# get cell model passport
file_CMP_ann <- "datasets/CMP_PORTAL/model_annotation/model_list_20230923.csv"
CMP_table <- read_csv(file_CMP_ann) %>%
  select(model_id, sample_id, model_name, synonyms, tissue, cancer_type, 
        tissue_status, COSMIC_ID, BROAD_ID, CCLE_ID) %>%
  rename(model_id_CMP = model_id, 
         sample_id_CMP = sample_id, 
         model_name_CMP = model_name)

CL_data_meta_QC <- left_join(by = "model_id_CMP", 
                             CL_data_meta_QC, 
                             CMP_table) %>%
  mutate(model_name_CMP_library = paste(model_name_CMP,library, sep = "_"))


# add results from additional screenings
fold_newscreens <- "datasets/PROJECT_SCORE/CLs_FROM_SANGER_09082023/"
files <- dir(fold_newscreens, 
             recursive = TRUE, 
             full.names = FALSE, 
             pattern = "\\count_matrix.tsv$")
CL_new_name <- str_split_i(files, pattern = "/", i = 1) 
CL_new_name <- setdiff(CL_new_name, "HT-29") # do not include HT-29

CL_data_new <- data.frame(
  location = NULL,
  complete_name = NULL,
  model_name_file = NULL,
  rep = NULL,
  model_name = NULL,
  passed_low_QC = NULL, 
  passed_high_QC = NULL,
  library = NULL
)
for (i in 1:length(CL_new_name)) {
  
  CL_tmp <- CL_new_name[i]
  files_tmp <- files[grepl(CL_tmp, files)]
  
  colnames_table <- colnames(read_table(paste0(fold_newscreens, files_tmp), 
                                        show_col_types = FALSE))
  rep_tmp <- colnames_table[!colnames_table %in% c("sgRNA", "gene", "HumanCRISPR_V1.1Plas")]
  complete_name_tmp <- paste(CL_tmp, rep_tmp, sep = "_")
  
  CL_data_new <- rbind(CL_data_new, data.frame(
    location = files_tmp,
    complete_name = complete_name_tmp,
    model_name_file = CL_tmp, 
    rep = rep_tmp,
    model_name = CL_tmp, 
    passed_low_QC = NA, 
    passed_high_QC = NA,
    library = "V1.1"
  ))
}
CL_data_new$model_name_CMP <- CL_data_new$model_name
CL_data_new_meta <- left_join(by = "model_name_CMP", 
                              CL_data_new, 
                              CMP_table) %>%
  mutate(model_name_CMP_library = paste(model_name_CMP,library, sep = "_"), 
         origin_data = "SANGER_SCREEN_09082023")

CL_data_meta_tot <- bind_rows(CL_data_meta_QC, CL_data_new_meta)

# save complete list
write_tsv(x = CL_data_meta_tot,
          file = sprintf("%sSINGLE_CL/model_list.tsv", fold), 
          col_names = T)

### save single files
CL_names <- unique(CL_data_meta_tot$model_name_CMP_library)
plasmid_names <- c("CRISPR_C6596666.sample", "ERS717283.plasmid")

for (CL in CL_names) {
  
  print(CL)
  CL_meta_curr <- CL_data_meta_tot %>% filter(model_name_CMP_library %in% CL)
  
  if (any(CL_meta_curr$origin_data != "PROJECT_SCORE_RELEASE_1")) {
    
    print("Already combined, copy table")
    CL_curr <- read_table(paste0(fold_newscreens, CL_meta_curr$location[1]), 
                      show_col_types = FALSE) 
    write.table(
      x = CL_curr, 
      file = sprintf("%sSINGLE_CL/%s_counts.tsv", fold, CL), 
      sep = "\t", 
      col.names = T, 
      row.names = F, 
      quote = F
    )
  }else{
    
    n_rep <- nrow(CL_meta_curr)
    tmp <- list()
    for (i in seq_len(n_rep)) {
      print(i)
      tmp[[i]] <- read_table(gzfile(paste0(fold, CL_meta_curr$location[i])), 
                             show_col_types = FALSE)
    }
    
    all_columns <- purrr::reduce(lapply(tmp, names), union)
    curr_plasmid <- intersect(plasmid_names, all_columns)
    common_columns <- c("sgRNA", "gene", curr_plasmid)
    CL_curr <- purrr::reduce(tmp, dplyr::inner_join, by = common_columns) %>%
      relocate(!!(common_columns)) # first 3 col are sgRNA, gene, and control
    
    print(colnames(CL_curr))
    
    write.table(
      x = CL_curr, 
      file = sprintf("%sSINGLE_CL/%s_counts.tsv", fold, CL), 
      sep = "\t", 
      col.names = T, 
      row.names = F, 
      quote = F
    )
  }
}

