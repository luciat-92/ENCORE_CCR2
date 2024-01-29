# validate batch correction over original values
library(renv)
renv::activate()
print(.libPaths())

# Combine all results across cell lines
library(CRISPRcleanRatSquared)
library(tidyverse)
library(ggpubr)
library(pheatmap)
library(reshape2)
library(R.utils)
library(argparse)
pdf(NULL) # avoid Rplots.pdf from being generated

parser <- ArgumentParser(description = "Summarise results from CRISPRcleanR^2 for all CLs")
parser$add_argument("--fold_original", type = "character", help = "path to folder output for no preprocessing")
parser$add_argument("--fold_batchcorr", type = "character", help = "path to folder output for batch corrected version")
parser$add_argument("--fold_output", type = "character", help = "path to folder to save output")

args <- parser$parse_args()
fold_original <- args$fold_original
fold_batchcorr <- args$fold_batchcorr
fold_output <- args$fold_output

########
# root_path <- "/group/iorio/lucia/"
# fold_output <- sprintf("%sCRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/SUMMARY_ORIGINAL_BATCH_CORRECTED/cavg/", root_path)
# fold_original <- sprintf("%sCRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/ORIGINAL/cavg/ALL_CLs/", root_path)
# fold_batchcorr <- sprintf("%sCRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/BATCH_CORRECTED/cavg/ALL_CLs/", root_path)
########

#################
### FUNCTIONS ###
#################
# compare
get_batch <- function(ID_lib){
  lib <- sort(sapply(str_split_1(ID_lib, pattern = ","), function(x) str_split_i(x, pattern = "_", 2)))
  lib <- paste0(lib, collapse = ",")
  return(lib)
}

get_match_correction <- function(original, batchcorr){
  
  batchcorr$batch <- unname(sapply(batchcorr$ID_lib, get_batch))
  batchcorr <- batchcorr %>%
    dplyr::mutate(ID_lib_CL = paste(ID_lib, CL, sep = "_")) %>%
    dplyr::select(ID_lib_CL, batch, info, info_subtype, Gene_Pair, sgRNA_ID_pair, avgFC, correction, correctedFC)
  
  original$batch <- unname(sapply(original$ID_lib, get_batch))
  original <- original %>%
    dplyr::mutate(ID_lib_CL = paste(ID_lib, CL, sep = "_")) %>%
    dplyr::select(ID_lib_CL, batch, info, info_subtype, Gene_Pair, sgRNA_ID_pair, avgFC, correction, correctedFC)
  
  tot <- dplyr::inner_join(original, batchcorr, 
                           by = c("ID_lib_CL", "batch","info", "info_subtype", "Gene_Pair", "sgRNA_ID_pair"), 
                           suffix = c("_original", "_batchcorr"))
  
  return(tot)
}
#################

file_original <- sprintf("%slogFC_sgRNA_CCR2correction.txt.gz", fold_original)
file_batchcorr <- sprintf("%slogFC_sgRNA_CCR2correction.txt.gz", fold_batchcorr)

out_original_all <- readr::read_tsv(gzfile(file_original), 
                                col_types = readr::cols(.default = "?",
                                                        sgRNA1_WGE_ID = "c", 
                                                        sgRNA2_WGE_ID = "c"))
out_batchcorr_all <- readr::read_tsv(gzfile(file_batchcorr), 
                                 col_types = readr::cols(.default = "?",
                                                         sgRNA1_WGE_ID = "c", 
                                                         sgRNA2_WGE_ID = "c"))

# check differences in correction distribution across all samples
names_CL_lib <- unique(out_original_all$CL_lib)
kruskal_results <- list()
print("KW test per CL")
for (i in 1:length(names_CL_lib)) {
  
  print(i)
  id <- names_CL_lib[i]
  # subset
  out_original <- out_original_all[out_original_all$CL_lib == id,]
  out_batchcorr <- out_batchcorr_all[out_batchcorr_all$CL_lib == id,]
  df_full <- get_match_correction(out_original, out_batchcorr)
  
  df_boxplot <- data.frame(correction = c(df_full$correction_original, 
                                          df_full$correction_batchcorr), 
                           batch = rep(df_full$batch, 2), 
                           type = c(rep("No Preprocessing", nrow(df_full)), 
                                    rep("ComBat Corrected", nrow(df_full)))) %>%
    dplyr::filter(!is.na(correction))
  
  if (out_original$lib[1] == "COLO") {
    keep_batch <- c("colo1", "colo2", "colo3")
  }
  if (out_original$lib[1] == "BRCA") {
    keep_batch <- c("brca1", "brca2", "brca3")
  }
  
  df_stat <- df_boxplot %>% 
    dplyr::filter(batch %in% keep_batch) %>%
    dplyr::group_by(batch, type) %>%
    dplyr::summarise(lim_low = quantile(correction, 0.25)  - 1.5*stats::IQR(correction), 
                     lim_up = quantile(correction, 0.75)  + 1.5*stats::IQR(correction))
  
  # Compute the Kruskal-Wallis test for each combination of type and batch
  kruskal_results[[i]] <- df_boxplot %>%
    dplyr::filter(batch %in% keep_batch) %>%
    dplyr::group_by(type) %>%
    dplyr::summarize(p_value = kruskal.test(x = correction, g = batch)$p.value, 
                     y_pos = quantile(correction, 0.75)  + 1.5*stats::IQR(correction) + 0.015,  
                     x_pos = keep_batch[grepl("2", keep_batch)]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(CL_lib = id)
  
}
kruskal_results <- do.call(rbind, kruskal_results)
# save 
write.table(x = kruskal_results %>% select(-y_pos, -x_pos), 
            file = sprintf("%scorrection_VS_batch_KWtest_perCL.txt", fold_output), 
            quote = F, 
            sep = "\t")
kruskal_results <- kruskal_results %>% 
  dplyr::mutate(p_value = case_when(p_value == 0 ~ 1e-300, 
                          TRUE ~ p_value)) %>%
  dplyr::mutate(log10p = -log10(p_value))

df <- dcast(data = kruskal_results, formula = CL_lib ~ type, value.var = "log10p") %>%
  dplyr::mutate(lib = str_split_i(CL_lib, "_", 2))

# save plot
pl <- ggplot(df, 
       aes(y = get("ComBat Corrected"), x = get("No Preprocessing"), color = lib)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("No Preprocessing") + 
  ylab("ComBat Corrected") + 
  ggtitle("-log10(P) Kruskal-Wallis test: CCR^2 correction ~ batch")

ggsave(filename = sprintf("%scorrection_VS_batch_KWtest_perCL.pdf", fold_output), 
       plot = pl, width = 5.2, height = 5)


# run across all CLs
print("KW test per tissue")
names_lib <- unique(out_original_all$lib)
kruskal_results <- list()
for (i in 1:length(names_lib)) {
  
  print(i)
  id <- names_lib[i]
  out_original <- out_original_all[out_original_all$lib == id,]
  out_batchcorr <- out_batchcorr_all[out_batchcorr_all$lib == id,]
  df_full <- get_match_correction(out_original, out_batchcorr)
  
  df_boxplot <- data.frame(correction = c(df_full$correction_original, 
                                          df_full$correction_batchcorr), 
                           batch = rep(df_full$batch, 2), 
                           type = c(rep("No Preprocessing", nrow(df_full)), 
                                    rep("ComBat Corrected", nrow(df_full)))) %>%
    dplyr::filter(!is.na(correction))
  
  if (out_original$lib[1] == "COLO") {
    keep_batch <- c("colo1", "colo2", "colo3")
  }
  if (out_original$lib[1] == "BRCA") {
    keep_batch <- c("brca1", "brca2", "brca3")
  }
  
  df_stat <- df_boxplot %>% 
    dplyr::filter(batch %in% keep_batch) %>%
    dplyr::group_by(batch, type) %>%
    dplyr::summarise(lim_low = quantile(correction, 0.25)  - 1.5*stats::IQR(correction), 
                     lim_up = quantile(correction, 0.75)  + 1.5*stats::IQR(correction))
  
  # Compute the Kruskal-Wallis test for each combination of type and batch
  kruskal_results[[i]] <- df_boxplot %>%
    dplyr::filter(batch %in% keep_batch) %>%
    dplyr::group_by(type) %>%
    dplyr::summarize(p_value = kruskal.test(x = correction, g = batch)$p.value, 
                     y_pos = quantile(correction, 0.75)  + 1.5*stats::IQR(correction) + 0.015,  
                     x_pos = keep_batch[grepl("2", keep_batch)]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(CL_lib = id)
}
kruskal_results <- do.call(rbind, kruskal_results)
# save 
write.table(x = kruskal_results, 
            file = sprintf("%scorrection_VS_batch_KWtest_perlib.txt", fold_output), 
            quote = F, 
            sep = "\t")

# plot model performances
print("Differences in model performances")

file_original <- sprintf("%sALL_CLs_avg_model_perf_DoubleGenePairs.tsv", fold_original)
file_batchcorr <- sprintf("%sALL_CLs_avg_model_perf_DoubleGenePairs.tsv", fold_batchcorr)

model_original <- readr::read_tsv(file_original)
model_batchcorr <- readr::read_tsv(file_batchcorr)

model_tot <- dplyr::inner_join(model_original, model_batchcorr, 
                               by = c("CL", "lib"),
                               suffix = c("_original", "_batchcorr"))
# annotate
pl <- ggplot(model_tot, 
             aes(x = mean_corr_original, y = mean_corr_batchcorr, color = lib)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  annotate("text", x = -Inf, y = Inf, 
           label = sprintf("ComBat Corr >= no preproc: %s", sum(model_tot$mean_corr_batchcorr >= model_tot$mean_corr_original)), 
           hjust = -0.1, vjust = 2, size = 4) +
  annotate("text", x = Inf, y = -Inf, 
           label = sprintf("ComBat Corr < no preproc: %s", sum(model_tot$mean_corr_batchcorr < model_tot$mean_corr_original)), 
           hjust = 1.1, vjust = -1, size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("No Preprocessing") + 
  ylab("ComBat Corrected") + 
  ggtitle("Pears. Corr. Single VS Pseudo single") 
pl
ggsave(filename = sprintf("%smodel_perf.pdf", fold_output), 
       plot = pl, width = 5.2, height = 5)

# plot gene selection
print("Differences in gene selection")
file_original <- sprintf("%sselection_gene_letANDsyn.txt", fold_original)
file_batchcorr <- sprintf("%sselection_gene_letANDsyn.txt", fold_batchcorr)

sel_original <- readr::read_tsv(file_original) %>%
  dplyr::mutate(preproc = "No Preprocessing")

sel_batchcorr <- readr::read_tsv(file_batchcorr) %>%
  dplyr::mutate(preproc = "ComBat Corrected")

# tot
sel_tot <- bind_rows(sel_original, sel_batchcorr) %>%
  dplyr::mutate(type_new = case_when(common == "in common" ~ "in common", 
                                     common == "unique" & type == "Uncorrected" ~ "unique (Uncorrected)", 
                                     common == "unique" & type == "Corrected" ~ "unique (CCR^2 Corrected)")) %>%
  dplyr::mutate(lib = factor(lib, levels = c("COLO", "BRCA")))

pl <- ggplot(sel_tot, 
             aes(x = type_new, fill = preproc)) +
  geom_bar(position = position_dodge()) + 
  facet_wrap(.~lib, ncol = 1) + 
  theme_bw() + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  coord_flip() + 
  xlab("") + 
  ylab("N. selected gene pairs across all CLs") 
pl
ggsave(filename = sprintf("%scount_selected_gene_pairs.pdf", fold_output), 
       plot = pl, width = 5.2, height = 3)

pl <- ggplot(sel_tot, 
             aes(x = CL, fill = preproc)) +
  geom_bar(position = position_dodge()) + 
  facet_wrap(.~lib, ncol = 1, scales = "free_y") + 
  theme_bw() + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  coord_flip() + 
  xlab("") + 
  ylab("N. selected gene pairs across all CLs") 
pl
ggsave(filename = sprintf("%scount_selected_gene_pairs_perCLs.pdf", fold_output), 
       plot = pl, width = 5.5, height = 5.5)




