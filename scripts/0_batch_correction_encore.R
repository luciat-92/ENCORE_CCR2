library(tidyverse)
library(sva)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(pROC)
library(CRISPRcleanR)
library(effectsize)
# setwd("~/workspace/ENCORE_CCR2/")
source("R/0_batchcorr_auxilary_functions.R")

#### 
setwd("/group/iorio/lucia/")
# set up folder
fold_data <- "datasets/ENCORE_SAMPLES_COPYNUMBER/"
fold_cmp <- "datasets/CMP_PORTAL/"
fold_input <- sprintf("%sDATA_FREEZE_v4_NOV_2023/cas9neg_rdsFiles/", fold_data)
fold_lib <- sprintf("%sDATA_META/LIBS/", fold_data)
# create output folder
system("mkdir -p CRISPR_combinatorial/data/encore/DATA_FREEZE_v4_NOV_2023/BATCH_CORRECTED/")
system("mkdir -p CRISPR_combinatorial/data/encore/DATA_FREEZE_v4_NOV_2023/ORIGINAL/")
fold_out <- "CRISPR_combinatorial/data/encore/DATA_FREEZE_v4_NOV_2023/"
####

lib_name <- c("COLO1", "COLO2", "COLO3", 
              "BRCA1", "BRCA2", "BRCA3") %>% tolower()

cneg_label <- c("c91", "c92")

res_cneg <- lapply(cneg_label, function(x) 
  load_data_rds(
    fold_input = fold_input, 
    fold_lib = fold_lib, 
    negcontrol = x)
  )

# get avg between the two plasmids
res_cneg_avg <- get_combined_cneg(
  res_cneg = res_cneg
)
names(res_cneg) <- cneg_label
res_cneg <- purrr::list_modify(res_cneg, cavg = res_cneg_avg)
cneg_label <- names(res_cneg)

###############################
##### NO BATCH CORRECTION #####
###############################

for (idx in 1:length(res_cneg)) {
  
  print(cneg_label[idx])
  
  model_encore_table <- get_sample_ann(
    data = res_cneg[[idx]]$data, 
    file_cmp = sprintf("%smodel_annotation/model_list_20230801.csv", fold_cmp))
  
  # get matrix
  data_COLO <- harmonize_per_CL(
    list_df = res_cneg[[idx]]$data[grepl("colo",names(res_cneg[[idx]]$data))], 
    CL_ann = model_encore_table)
  
  data_BRCA <- harmonize_per_CL(
    list_df = res_cneg[[idx]]$data[grepl("brca", names(res_cneg[[idx]]$data))], 
    CL_ann = model_encore_table)
  
  # get final table
  df_COLO <- get_complete_table(
    list_df = res_cneg[[idx]]$data[grepl("colo",names(res_cneg[[idx]]$data))], 
    list_matrix = data_COLO)
  
  df_BRCA <- get_complete_table(
    list_df = res_cneg[[idx]]$data[grepl("brca",names(res_cneg[[idx]]$data))], 
    list_matrix = data_BRCA)
  
  # save output
  write.table(file = sprintf("%sORIGINAL/COLO_%s_FINAL_EXACT_logFC_sgRNA.txt", fold_out, cneg_label[idx]), 
              x = df_COLO, 
              quote = F, 
              col.names = T, 
              row.names = F, 
              sep = "\t")
  
  write.table(file = sprintf("%sORIGINAL/BRCA_%s_FINAL_EXACT_logFC_sgRNA.txt", fold_out, cneg_label[idx]), 
              x = df_BRCA, 
              quote = F, 
              col.names = T, 
              row.names = F, 
              sep = "\t")
  
  # save combined libraries (only once!)
  if (!file.exists(sprintf("%sORIGINAL/ENCORE_GI_COREAD_Library_ALL.txt", fold_out))) {
    write.table(file = sprintf("%sORIGINAL/ENCORE_GI_COREAD_Library_ALL.txt", fold_out), 
                x = bind_rows(res_cneg[[1]]$library[grepl("colo",names(res_cneg[[1]]$library))]), 
                quote = F, 
                col.names = T, 
                row.names = F, 
                sep = "\t")
  }
  
  if (!file.exists(sprintf("%sORIGINAL/ENCORE_GI_BRCA_Library_ALL.txt", fold_out))) {
    write.table(file = sprintf("%sORIGINAL/ENCORE_GI_BRCA_Library_ALL.txt", fold_out), 
              x = bind_rows(res_cneg[[1]]$library[grepl("brca",names(res_cneg[[1]]$library))]), 
              quote = F, 
              col.names = T, 
              row.names = F, 
              sep = "\t")
  }
}

############################
##### BATCH CORRECTION #####
############################
plots_batch_correction <- TRUE


# for (idx in 1:length(res_cneg)) {
  idx=1
  id_colo <- grepl("colo",names(res_cneg[[idx]]$data))
  id_brca <- grepl("brca",names(res_cneg[[idx]]$data))
  
  model_encore_table <- get_sample_ann(
    data = res_cneg[[idx]]$data, 
    file_cmp = sprintf("%smodel_annotation/model_list_20230801.csv", fold_cmp))
  
  # plots
  if (plots_batch_correction) {
  
    common_COLO <- pca_commonpairs_function(
      list_df = res_cneg[[idx]]$data[id_colo], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/COLO_%s_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      show_plot = TRUE, 
      save_ext = "pdf"
    ) 
    
    common_BRCA <- pca_commonpairs_function(
      list_df = res_cneg[[idx]]$data[id_brca], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/BRCA_%s_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      show_plot = TRUE, 
      save_ext = "pdf"
    )
    
    plot_dist_commonpairs(
      list_df = res_cneg[[idx]]$data[id_colo], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/COLO_%s_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      save_ext = "pdf"
    )
    
    plot_dist_commonpairs(
      list_df = res_cneg[[idx]]$data[id_brca], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/BRCA_%s_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      save_ext = "pdf"
    )
    
    plot_dist_PPV(
      list_df = res_cneg[[idx]]$data[id_colo],
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/COLO_%s_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      save_ext = "pdf"
    )

    plot_dist_PPV(
      list_df = res_cneg[[idx]]$data[id_brca], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/BRCA_%s_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      save_ext = "pdf"
    )
    
  }
  
  
  
#}












