library(tidyverse)
library(sva)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(pROC)
library(CRISPRcleanR)
library(effectsize)

source("R/0_auxilary_functions.R")

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
  system(paste0("mkdir -p ", fold_out, "/ORIGINAL/", cneg_label[idx]))
  write.table(file = sprintf("%sORIGINAL/%s/COLO_FINAL_EXACT_logFC_sgRNA.txt", fold_out, cneg_label[idx]),
              x = df_COLO,
              quote = F,
              col.names = T,
              row.names = F,
              sep = "\t")

  write.table(file = sprintf("%sORIGINAL/%s/BRCA_FINAL_EXACT_logFC_sgRNA.txt", fold_out, cneg_label[idx]),
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
data("ADaM2021_essential")
plots_batch_correction <- TRUE

for (idx in 1:length(res_cneg)) {
  
  print(paste0("########## ", cneg_label[idx], " ###########"))
  id_colo <- grepl("colo",names(res_cneg[[idx]]$data))
  id_brca <- grepl("brca",names(res_cneg[[idx]]$data))
  
  model_encore_table <- get_sample_ann(
    data = res_cneg[[idx]]$data, 
    file_cmp = sprintf("%smodel_annotation/model_list_20230801.csv", fold_cmp))
  
  system(paste0("mkdir -p ", fold_out, "/BATCH_CORRECTED/", cneg_label[idx]))
  # plots
  if (plots_batch_correction) {
  
    common_COLO <- pca_commonpairs_function(
      list_df = res_cneg[[idx]]$data[id_colo], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/%s/COLO_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      show_plot = TRUE, 
      save_ext = "pdf"
    ) 
    
    common_BRCA <- pca_commonpairs_function(
      list_df = res_cneg[[idx]]$data[id_brca], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/%s/BRCA_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      show_plot = TRUE, 
      save_ext = "pdf"
    )
    
    plot_dist_commonpairs(
      list_df = res_cneg[[idx]]$data[id_colo], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/%s/COLO_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      save_ext = "pdf"
    )
    
    plot_dist_commonpairs(
      list_df = res_cneg[[idx]]$data[id_brca], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/%s/BRCA_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      save_ext = "pdf"
    )
    
    plot_dist_PPV(
      list_df = res_cneg[[idx]]$data[id_colo],
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/%s/COLO_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      save_ext = "pdf"
    )

    plot_dist_PPV(
      list_df = res_cneg[[idx]]$data[id_brca], 
      CL_ann = model_encore_table,
      outfold = sprintf("%sBATCH_CORRECTED/%s/BRCA_", fold_out, cneg_label[idx]), 
      save_plot = TRUE, 
      save_ext = "pdf"
    )
    
  }
  
  # which kNN should be used?
  validation_COLO <- validate_NN_approximation(
    list_df = res_cneg[[idx]]$data[id_colo], 
    CL_ann = model_encore_table,
    outfold = sprintf("%sBATCH_CORRECTED/%s/COLO_", fold_out, cneg_label[idx]), 
    save_plot = TRUE, 
    save_ext = "pdf"
  )
  
  validation_BRCA <- validate_NN_approximation(
    list_df = res_cneg[[idx]]$data[id_brca], 
    CL_ann = model_encore_table,
    outfold = sprintf("%sBATCH_CORRECTED/%s/BRCA_", fold_out, cneg_label[idx]), 
    save_plot = TRUE, 
    save_ext = "pdf"
  )
  
  best_kNN_COLO <- validation_COLO$cohensd_kNN %>% 
    group_by(kNN) %>% 
    summarise(median_D = median(cohens_d)) %>%
    filter(median_D == min(median_D)) %>%
    pull(kNN)
  print(paste0("best kNN for COLO based on cohens D: ", best_kNN_COLO))
  
  best_kNN_BRCA <- validation_BRCA$cohensd_kNN %>% 
    group_by(kNN) %>% 
    summarise(median_D = median(cohens_d)) %>%
    filter(median_D == min(median_D)) %>%
    pull(kNN)
  print(paste0("best kNN for BRCA based on cohens D: ", best_kNN_BRCA))
  
  # get all corrected dataset
  data_COLO <- adjust_alldata_kNN(
    list_df = res_cneg[[idx]]$data[id_colo], 
    CL_ann = model_encore_table,
    kNN = best_kNN_COLO, 
    outfold = sprintf("%sBATCH_CORRECTED/%s/COLO_", fold_out, cneg_label[idx]), 
    save_plot = TRUE, 
    show_plot = TRUE, 
    save_ext = "pdf") 
  
  # get all corrected dataset
  data_BRCA <- adjust_alldata_kNN(
    list_df = res_cneg[[idx]]$data[id_brca], 
    CL_ann = model_encore_table,
    kNN = best_kNN_BRCA, 
    outfold = sprintf("%sBATCH_CORRECTED/%s/BRCA_", fold_out, cneg_label[idx]), 
    save_plot = TRUE, 
    show_plot = TRUE, 
    save_ext = "pdf") 
  
  # plot distribution
  COLO_allCLs <- plot_CL_distribution(original = data_COLO$original, 
                                      adjusted = data_COLO$adj, 
                                      list_df = res_cneg[[idx]]$data[id_colo],
                                      common_pairs = data_COLO$combat$common_pairs, 
                                      outfold = sprintf("%sBATCH_CORRECTED/%s/COLO_", fold_out, cneg_label[idx]), 
                                      save_plot = TRUE, 
                                      show_plot = TRUE, 
                                      save_ext = "png") 
  
  BRCA_allCLs <- plot_CL_distribution(original = data_BRCA$original, 
                                      adjusted = data_BRCA$adj, 
                                      list_df = res_cneg[[idx]]$data[id_brca],
                                      common_pairs = data_BRCA$combat$common_pairs, 
                                      outfold = sprintf("%sBATCH_CORRECTED/%s/BRCA_", fold_out, cneg_label[idx]), 
                                      save_plot = TRUE, 
                                      show_plot = TRUE,
                                      save_ext = "png")
  
  # get final table
  df_adj_COLO <- get_complete_table(
    list_df = res_cneg[[idx]]$data[grepl("colo",names(res_cneg[[idx]]$data))], 
    list_matrix = data_COLO$adj)
  
  df_or_COLO <- get_complete_table(
    list_df = res_cneg[[idx]]$data[grepl("colo",names(res_cneg[[idx]]$data))], 
    list_matrix = data_COLO$original)
  
  df_adj_BRCA <- get_complete_table(
    list_df = res_cneg[[idx]]$data[grepl("brca",names(res_cneg[[idx]]$data))], 
    list_matrix = data_BRCA$adj)
  
  df_or_BRCA <- get_complete_table(
    list_df = res_cneg[[idx]]$data[grepl("brca",names(res_cneg[[idx]]$data))], 
    list_matrix = data_BRCA$original)

  # save output
  write.table(file = sprintf("%sBATCH_CORRECTED/%s/COLO_FINAL_EXACT_logFC_sgRNA.txt", fold_out, cneg_label[idx]), 
              x = df_adj_COLO, 
              quote = F, 
              col.names = T, 
              row.names = F, 
              sep = "\t")
  
  write.table(file = sprintf("%sBATCH_CORRECTED/%s/BRCA_FINAL_EXACT_logFC_sgRNA.txt", fold_out, cneg_label[idx]), 
              x = df_adj_BRCA, 
              quote = F, 
              col.names = T, 
              row.names = F, 
              sep = "\t")
  
  # save combined libraries (only once!)
  if (!file.exists(sprintf("%sBATCH_CORRECTED/ENCORE_GI_COREAD_Library_ALL.txt", fold_out))) {
    write.table(file = sprintf("%sBATCH_CORRECTED/ENCORE_GI_COREAD_Library_ALL.txt", fold_out), 
                x = bind_rows(res_cneg[[1]]$library[grepl("colo",names(res_cneg[[1]]$library))]), 
                quote = F, 
                col.names = T, 
                row.names = F, 
                sep = "\t")
  }
  
  if (!file.exists(sprintf("%sBATCH_CORRECTED/ENCORE_GI_BRCA_Library_ALL.txt", fold_out))) {
    write.table(file = sprintf("%sBATCH_CORRECTED/ENCORE_GI_BRCA_Library_ALL.txt", fold_out), 
                x = bind_rows(res_cneg[[1]]$library[grepl("brca",names(res_cneg[[1]]$library))]), 
                quote = F, 
                col.names = T, 
                row.names = F, 
                sep = "\t")
  }
  
  # save parameters for each guide pair
  param <- data_COLO$param_all
  save(param, 
       file = sprintf("%sBATCH_CORRECTED/%s/COLO_FINAL_EXACT_logFC_sgRNA_ComBatParam.RData", fold_out, cneg_label[idx]))
  
  param <- data_BRCA$param_all
  save(param, 
       file = sprintf("%sBATCH_CORRECTED/%s/BRCA_FINAL_EXACT_logFC_sgRNA_ComBatParam.RData", fold_out, cneg_label[idx]))
  
  
  ## external validation: check CL specific distributions before and after correction ###
  ktest_COLO <- test_distributions_per_class(data_adj = df_adj_COLO, 
                                             data_or = df_or_COLO, 
                                             outfold = sprintf("%sBATCH_CORRECTED/%s/COLO_", fold_out, cneg_label[idx]), 
                                             save_plot = TRUE, 
                                             show_plot = TRUE, 
                                             save_ext = "pdf")
  ktest_COLO_all <- test_distributions_per_class_allCLs(COLO_allCLs, 
                                      outfold = sprintf("%sBATCH_CORRECTED/%s/COLO_", fold_out, cneg_label[idx]), 
                                      save_plot = TRUE, 
                                      save_ext = "pdf")
  
  ktest_BRCA <- test_distributions_per_class(data_adj = df_adj_BRCA, 
                                             data_or = df_or_BRCA, 
                                             outfold = sprintf("%sBATCH_CORRECTED/%s/BRCA_", fold_out, cneg_label[idx]), 
                                             save_plot = TRUE, 
                                             show_plot = TRUE, 
                                             save_ext = "pdf")
  ktest_BRCA_all <- test_distributions_per_class_allCLs(BRCA_allCLs, 
                                                        outfold = sprintf("%sBATCH_CORRECTED/%s/BRCA_", fold_out, cneg_label[idx]), 
                                                        save_plot = TRUE, 
                                                        save_ext = "pdf")
  
  # those in library singletons, do they have an in balance in essential genes?
  COLO_lib_genes <- plot_library_genes(data_adj = df_adj_COLO, 
                                       data_or = df_or_COLO, 
                                       essential_genes = ADaM2021_essential,  
                                       outfold = sprintf("%sBATCH_CORRECTED/%s/COLO_", fold_out, cneg_label[idx]), 
                                       save_plot = TRUE, 
                                       show_plot = TRUE, 
                                       save_ext = "pdf")
  
  BRCA_lib_genes <- plot_library_genes(data_adj = df_adj_BRCA, 
                                       data_or = df_or_BRCA, 
                                       essential_genes = ADaM2021_essential,  
                                       outfold = sprintf("%sBATCH_CORRECTED/%s/BRCA_", fold_out, cneg_label[idx]), 
                                       save_plot = TRUE, 
                                       show_plot = TRUE, 
                                       save_ext = "pdf")
  
}












