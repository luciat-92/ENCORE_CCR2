# Apply systematically across all CLs
library(renv)
renv::activate()
print(.libPaths())

# Combine all results across cell lines
library(CRISPRcleanRatSquared)
library(CRISPRcleanR)
library(S4Vectors)
library(tidyverse)
library(ggpubr)
library(pheatmap)
library(reshape2)
library(R.utils)
library(argparse)
pdf(NULL) # avoid Rplots.pdf from being generated

parser <- ArgumentParser(description = "Summarise results from CRISPRcleanR^2 for all CLs")
parser$add_argument("--fold_input_dual", type = "character", help = "path to folder with preprocessed data from ENCORE")
parser$add_argument("--fold_output", type = "character", help = "path to folder to save output")
parser$add_argument("--fold_CN", type = "character", help = "path to folder with CN file")
parser$add_argument("--root_path", type = "character", help = "base path, default is group-share", default = "/group/iorio/lucia/")

args <- parser$parse_args()
fold_input_dual <- args$fold_input_dual
fold_output <- args$fold_output
fold_CN <- args$fold_CN 
root_path <- args$root_path

source("R/3_auxilary_functions.R")

###########################################
# root_path <- "/group/iorio/lucia/"
# fold_input_dual <- sprintf("%sCRISPR_combinatorial/data/encore/DATA_FREEZE_v4_NOV_2023/ORIGINAL/c91/", root_path)
# fold_output <- sprintf("%sCRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/ORIGINAL/c91/", root_path)
# fold_CN <- sprintf("%sdatasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4_NOV_2023/METADATA_FEB2023/COPY_NUMBER/NEW_COPY_NUMBER/", root_path)
###########################################

##################################
copy_number_file <- sprintf("%sMERGED_SEGMENT_COPYNUMBER.txt", fold_CN)
system(sprintf("mkdir -p %s/ALL_CLs/",  fold_output))

model_encore_table <- read_tsv(sprintf("%smodel_list.tsv", fold_input_dual), 
                               show_col_types = FALSE) %>%
  dplyr::mutate(model_name_CMP_lib = paste(model_name_CMP, lib, sep = "_")) %>%
  dplyr::distinct(model_name_CMP_lib, .keep_all = TRUE) %>%
  dplyr::select(model_name_CMP_lib, model_name_CMP, lib, model_id_CMP)

output <- get_all_CLs(
  model_encore_table = model_encore_table, 
  fold = fold_output, 
  copy_number_file = copy_number_file)
# save
# save(output, file = sprintf("%s/ALL_CLs/complete_output.RData",  fold_output))
write.table(
  output$dual_FC, 
  file = sprintf("%s/ALL_CLs/logFC_sgRNA_CCR2correction.txt",  fold_output), 
  quote = F, 
  sep = "\t", 
  row.names = F, 
  col.names = T)
gzip(sprintf("%s/ALL_CLs/logFC_sgRNA_CCR2correction.txt", fold_output), 
     overwrite = TRUE) # compress

# load genome-wide single
single_gw <- get_all_CLs_singleGW(
  model_encore_table = model_encore_table, 
  fold = fold_output
)

gc();gc();gc()

####### plots #######
# Robject_ess <- "CRISPR_combinatorial/data/FiPer_outputs.RData" # from https://github.com/DepMap-Analytics/CoRe/blob/master/notebooks/data/preComputed/
# load(Robject_ess) # essential genes from FiPer method
# FiPer_essential <- Perc_AUC
data("ADaM2021_essential")
essential_genes <- ADaM2021_essential

##### plot model performances #####
## plot min thr and max thr
#thr_correction <- output$thr_correction
#plot(thr_correction$min_thr, thr_correction$min_correction, pch = 20)
#abline(a = 0, b = 1)
#plot(thr_correction$max_thr, thr_correction$max_correction, pch = 20)
#abline(a = 0, b = 1)

# plot model est p-values:
df_model_est <- output$model_est %>%
  dplyr::filter(type == "DOUBLE_CUT_PAIR") %>%
  dplyr::mutate(log10P = -log10(get("Pr(>|t|)"))) %>%
  dplyr::mutate(log10P = case_when(
    log10P == Inf ~ -log10(.Machine$double.xmin),
    TRUE ~ log10P))
rownames(df_model_est) <- NULL
pl <- ggplot(df_model_est, aes(x = CL_lib, fill = position, y = log10P)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  facet_wrap(.~feature, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("") +
  ylab("-log10(p-value)")
pl
ggsave(pl, filename = sprintf("%s/ALL_CLs/model_estimates_pvalue.png",  fold_output), width = 7, height = 5)


# split per CL
plot_model_perf(model_perf = output$model_perf,
                saveToFig = TRUE, 
                saveFormat = "png", 
                EXPname = "ALL_CLs", 
                outdir = sprintf("%s/ALL_CLs/",  fold_output))

# find relationship with model performance
df_model_perf <- plot_model_perf_vs_corr(
  dual_FC_withCNA = output$dual_FC_withCNA,
  single_FC_withCNA = output$single_FC_withCNA,
  single_FC_gw = single_gw$single_FC,
  model_perf = output$model_perf,
  tissue_name = "COLO",
  saveToFig = TRUE, 
  saveFormat = "png", 
  EXPname = "ALL_CLs", 
  outdir = sprintf("%s/ALL_CLs/", fold_output))

df_model_perf <- plot_model_perf_vs_corr(
  dual_FC_withCNA = output$dual_FC_withCNA,
  single_FC_withCNA = output$single_FC_withCNA,
  single_FC_gw = single_gw$single_FC,
  model_perf = output$model_perf,
  tissue_name = "BRCA",
  saveToFig = TRUE, 
  saveFormat = "png", 
  EXPname = "ALL_CLs", 
  outdir = sprintf("%s/ALL_CLs/", fold_output))

### summary ###
model_encore_table <- model_encore_table %>% arrange(lib)
CL <- unique(model_encore_table$model_name_CMP)
tab_count <- table(model_encore_table$model_name_CMP, model_encore_table$lib)
tab_count <- tab_count[match(CL, rownames(tab_count)),]
CL_summary <- pheatmap::pheatmap(tab_count,
                                 color = c("grey80", "black"),
                                 breaks = c(0, 0.5, 1),
                                 cluster_rows = TRUE,
                                 cluster_cols = FALSE,
                                 treeheight_row = 0)

save_pheatmap_png <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width = width, height = height, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(CL_summary,
                  sprintf("%s/ALL_CLs/CL_library_summary.png",  fold_output),
                  width = 2.5,
                  height = 5.5)

info_libraries <- output$dual_FC %>%
  group_by(ID_lib) %>%
  summarise(lib = unique(lib),
            info_subtype = unique(info_subtype),
            info = unique(info),
            Gene_Pair = unique(Gene_Pair),
            Gene1 =  unique(Gene1),
            Gene2 = unique(Gene2),
            sgRNA1_Library = unique(sgRNA1_Library),
            sgRNA2_Library = unique(sgRNA2_Library),
            sgRNA1_WGE_ID = unique(sgRNA1_WGE_ID),
            sgRNA2_WGE_ID = unique(sgRNA2_WGE_ID)) %>%
  ungroup() %>%
  mutate(tissue = lib)

pl <- ggplot(info_libraries, aes(x = lib, fill = tissue)) +
  geom_bar(stat = "count", width = 0.5) +
  scale_fill_manual(values = c("lightblue", "red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("") +
  ylab("N. of guide pairs")
pl
ggsave(pl, filename = sprintf("%s/ALL_CLs/count_guide_pairs_lib.png",  fold_output), width = 3, height = 3)

gene_info_libraries <- info_libraries %>%
  group_by(lib) %>%
  distinct(Gene_Pair, .keep_all = TRUE) %>%
  ungroup()

pl <- ggplot(gene_info_libraries, aes(x = lib, fill = tissue)) +
  geom_bar(stat = "count", width = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("lightblue", "red")) +
  xlab("") +
  ylab("N. of gene pairs")
pl
ggsave(pl, filename = sprintf("%s/ALL_CLs/count_gene_pairs_lib.png",  fold_output), width = 3, height = 3)

gene1_all_info_libraries <- info_libraries %>%
  group_by(tissue) %>%
  distinct(Gene1, .keep_all = TRUE) %>%
  rename(Gene = Gene1) %>%
  mutate(position = "Gene1")

gene2_all_info_libraries <- info_libraries %>%
  group_by(tissue) %>%
  distinct(Gene2, .keep_all = TRUE) %>%
  rename(Gene = Gene2) %>%
  mutate(position = "Gene2")
genepos_all_info_library <- rbind(gene1_all_info_libraries, gene2_all_info_libraries)

pl <- ggplot(genepos_all_info_library, aes(x = tissue, fill = tissue)) +
  geom_bar(stat = "count", width = 0.5) +
  facet_wrap(.~position, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("lightblue", "red")) +
  xlab("") +
  ylab("N. of unique genes")
pl
ggsave(pl, filename = sprintf("%s/ALL_CLs/count_genes_tissue.png",  fold_output), width = 3, height = 3)

####
# plot copy number from before correction logFC
plot_CN_vs_logFC(dual_FC_withCNA = output$dual_FC_withCNA,
                 saveToFig = TRUE,
                 saveFormat = "png",
                 EXPname = "ALL_CLs",
                 outdir = sprintf("%s/ALL_CLs/",  fold))

plot_CN_vs_logFC(dual_FC_withCNA = output$dual_FC_withCNA[grepl("COLO", output$dual_FC_withCNA$lib),],
                 saveToFig = TRUE,
                 saveFormat = "png",
                 EXPname = "ALL_CLs",
                 outdir = sprintf("%s/ALL_CLs/COLO_",  fold_output))

plot_CN_vs_logFC(dual_FC_withCNA = output$dual_FC_withCNA[grepl("BRCA", output$dual_FC_withCNA$lib),],
                 saveToFig = TRUE,
                 saveFormat = "png",
                 EXPname = "ALL_CLs",
                 outdir = sprintf("%s/ALL_CLs/BRCA_",  fold_output))

# plot number of pairs with Gene1 or Gene2 Amp per cell line
plot_count_amp(dual_FC_withCNA = output$dual_FC_withCNA,
               saveToFig = TRUE,
               saveFormat = "png",
               EXPname = "ALL_CLs",
               outdir = sprintf("%s/ALL_CLs/",  fold_output))

# count number of pairs
plotCNA_count(dual_FC_withCNA = output$dual_FC_withCNA,
              EXPname = "ALL_CLs",
              saveToFig = TRUE,
              saveFormat = "png",
              essential_genes = essential_genes,
              outdir = sprintf("%s/ALL_CLs/",  fold_output))

plotCNA_count(dual_FC_withCNA = output$dual_FC_withCNA[grepl("COLO", output$dual_FC_withCNA$lib),],
              EXPname = "ALL_CLs",
              saveToFig = TRUE,
              saveFormat = "png",
              essential_genes = essential_genes,
              outdir = sprintf("%s/ALL_CLs/COLO_",  fold_output))

plotCNA_count(dual_FC_withCNA = output$dual_FC_withCNA[grepl("BRCA", output$dual_FC_withCNA$lib),],
              EXPname = "ALL_CLs",
              saveToFig = TRUE,
              saveFormat = "png",
              essential_genes = essential_genes,
              outdir = sprintf("%s/ALL_CLs/BRCA_",  fold_output))

###### plot CN all #####
plotCNA(dual_FC_withCNA = output$dual_FC_withCNA, 
        EXPname = "ALL_CLs", 
        saveToFig = TRUE, 
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/",  fold_output))

plotCNA_diff(dual_FC_withCNA = output$dual_FC_withCNA, 
             EXPname = "ALL_CLs", 
             saveToFig = TRUE, 
             saveFormat = "png", 
             outdir = sprintf("%s/ALL_CLs/",  fold_output))

plotCNA(dual_FC_withCNA = output$dual_FC_withCNA, 
        EXPname = "ALL_CLs", 
        essential_genes = essential_genes,
        saveToFig = TRUE, 
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/",  fold_output))

plotCNA_diff(dual_FC_withCNA = output$dual_FC_withCNA, 
        EXPname = "ALL_CLs", 
        essential_genes = essential_genes,
        saveToFig = TRUE, 
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/",  fold_output))

plotCNA(dual_FC_withCNA = output$dual_FC_withCNA,
        EXPname = "ALL_CLs", 
        essential_genes = essential_genes,
        exclude_ess = FALSE, 
        saveToFig = TRUE, 
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/",  fold_output))

plotCNA_diff(dual_FC_withCNA = output$dual_FC_withCNA, 
             EXPname = "ALL_CLs", 
             essential_genes = essential_genes,
             exclude_ess = FALSE, 
             saveToFig = TRUE, 
             saveFormat = "png", 
             outdir = sprintf("%s/ALL_CLs/",  fold_output))

###### plot CN COLO #####
plotCNA(dual_FC_withCNA = output$dual_FC_withCNA[grepl("COLO", output$dual_FC_withCNA$lib),], 
        EXPname = "ALL_CLs", 
        saveToFig = TRUE, 
        essential_genes = essential_genes,
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/COLO_",  fold_output))

plotCNA_diff(dual_FC_withCNA = output$dual_FC_withCNA[grepl("COLO", output$dual_FC_withCNA$lib),], 
        EXPname = "ALL_CLs", 
        saveToFig = TRUE, 
        essential_genes = essential_genes,
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/COLO_",  fold_output))

plotCNA(dual_FC_withCNA = output$dual_FC_withCNA[grepl("COLO", output$dual_FC_withCNA$lib),], 
        EXPname = "ALL_CLs", 
        saveToFig = TRUE, 
        essential_genes = essential_genes,
        exclude_ess = FALSE, 
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/COLO_",  fold_output))

plotCNA_diff(dual_FC_withCNA = output$dual_FC_withCNA[grepl("COLO", output$dual_FC_withCNA$lib),], 
             EXPname = "ALL_CLs", 
             saveToFig = TRUE, 
             essential_genes = essential_genes,
             exclude_ess = FALSE, 
             saveFormat = "png", 
             outdir = sprintf("%s/ALL_CLs/COLO_",  fold_output))

plotCNA(dual_FC_withCNA = output$dual_FC_withCNA[grepl("COLO", output$dual_FC_withCNA$lib),], 
        EXPname = "ALL_CLs", 
        saveToFig = TRUE, 
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/COLO_",  fold_output))

plotCNA_diff(dual_FC_withCNA = output$dual_FC_withCNA[grepl("COLO", output$dual_FC_withCNA$lib),], 
             EXPname = "ALL_CLs", 
             saveToFig = TRUE, 
             saveFormat = "png", 
             outdir = sprintf("%s/ALL_CLs/COLO_",  fold_output))

###### plot CN BRCA #####
plotCNA(dual_FC_withCNA = output$dual_FC_withCNA[grepl("BRCA", output$dual_FC_withCNA$lib),], 
        EXPname = "ALL_CLs", 
        essential_genes = essential_genes,
        saveToFig = TRUE, 
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/BRCA_",  fold_output))

plotCNA_diff(dual_FC_withCNA = output$dual_FC_withCNA[grepl("BRCA", output$dual_FC_withCNA$lib),], 
             EXPname = "ALL_CLs", 
             essential_genes = essential_genes,
             saveToFig = TRUE, 
             saveFormat = "png", 
             outdir = sprintf("%s/ALL_CLs/BRCA_",  fold_output))

plotCNA(dual_FC_withCNA = output$dual_FC_withCNA[grepl("BRCA", output$dual_FC_withCNA$lib),], 
        EXPname = "ALL_CLs", 
        saveToFig = TRUE, 
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/BRCA_",  fold_output))

plotCNA_diff(dual_FC_withCNA = output$dual_FC_withCNA[grepl("BRCA", output$dual_FC_withCNA$lib),], 
             EXPname = "ALL_CLs", 
             saveToFig = TRUE, 
             saveFormat = "png", 
             outdir = sprintf("%s/ALL_CLs/BRCA_",  fold_output))


plotCNA(dual_FC_withCNA = output$dual_FC_withCNA[grepl("BRCA", output$dual_FC_withCNA$lib),], 
        EXPname = "ALL_CLs", 
        essential_genes = essential_genes,
        exclude_ess = FALSE, 
        saveToFig = TRUE, 
        saveFormat = "png", 
        outdir = sprintf("%s/ALL_CLs/BRCA_",  fold_output))

plotCNA_diff(dual_FC_withCNA = output$dual_FC_withCNA[grepl("BRCA", output$dual_FC_withCNA$lib),], 
             EXPname = "ALL_CLs", 
             essential_genes = essential_genes,
             exclude_ess = FALSE, 
             saveToFig = TRUE, 
             saveFormat = "png", 
             outdir = sprintf("%s/ALL_CLs/BRCA_",  fold_output))

##### selection of gene pairs ####
sel_gene <- output$selection_gene_letANDsyn %>%
  mutate(common_type = paste0(common, " (", type, ")")) %>%
  group_by(CL_lib) %>%
  distinct(Gene_Pair, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(common_type = case_when(
    grepl("unique", common_type) ~ common_type, 
    TRUE ~ "in common"))

count_sel_gene <- sel_gene %>% 
  group_by(CL_lib, common_type) %>%
  summarise(count_pairs = n(), 
            CL = unique(CL), 
            lib = unique(lib)) 

color_cl <- data.frame(tissue = c("COLO", "BRCA"), col = c("red", "lightblue"))
levels_colo <- count_sel_gene %>% 
  filter(grepl("COLO",lib)) %>%
  group_by(CL) %>% 
  summarise(tot = sum(count_pairs), 
            CL = unique(CL)) %>%
  arrange(desc(tot)) %>%
  pull(CL)

pl_colo <- ggplot(data = subset(count_sel_gene, grepl("COLO",lib)), aes(x = factor(CL, levels = levels_colo), 
                                  y = count_pairs,
                                  fill = common_type)) +
  geom_bar(stat = "identity", width = 0.8) + 
  facet_wrap(.~lib, nrow = 1) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank()) +
  ylab("N. selected gene pairs") + 
  coord_flip()

levels_brca <- count_sel_gene %>% 
  filter(grepl("BRCA",lib)) %>%
  group_by(CL) %>% 
  summarise(tot = sum(count_pairs), 
            CL = unique(CL)) %>%
  arrange(desc(tot)) %>%
  pull(CL)

pl_brca <- ggplot(data = subset(count_sel_gene, grepl("BRCA",lib)), aes(x = factor(CL, levels = levels_brca), 
                                                                        y = count_pairs,
                                                                        fill = common_type)) +
  geom_bar(stat = "identity", width = 0.8) + 
  facet_wrap(.~lib, nrow = 1) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank()) +
  ylab("N. selected gene pairs") + 
  coord_flip()
  
pl <- ggarrange(plotlist = list(pl_colo, pl_brca), ncol = 1, 
                common.legend = TRUE, 
                align = "v", 
                heights = c(1, 0.7))
print(pl)
ggsave(pl,  
       filename = sprintf("%s/ALL_CLs/ALL_CLs_count_selected_gene_pairs.png", fold_output), width = 5, height = 8)

# compare to amp count
count_amp_COLO <- output$dual_FC_withCNA %>%
  filter(grepl("COLO", lib)) %>%
  group_by(CL_lib) %>%
  summarise(n_amp = sum(as.character(Max_CN_cat) == "Amp"), 
            perc_amp = sum(as.character(Max_CN_cat) %in% c("Amp"))/length(Max_CN_cat), 
            mean_correction_dual = mean(correction, na.rm = TRUE))

count_novel_COLO <- count_sel_gene %>%
  filter(grepl("COLO",lib), grepl("unique", common_type)) %>%
  group_by(CL_lib) %>%
  summarise(n_novel = sum(count_pairs))

count_COLO <- full_join(count_amp_COLO, count_novel_COLO) %>%
  replace(is.na(.), 0) %>%
  mutate(tissue = "COLO")

count_amp_BRCA <- output$dual_FC_withCNA %>%
  filter(grepl("BRCA", lib)) %>%
  group_by(CL_lib) %>%
  summarise(n_amp = sum(as.character(Max_CN_cat) == "Amp"), 
            perc_amp = sum(as.character(Max_CN_cat) %in% c("Amp"))/length(Max_CN_cat), 
            mean_correction_dual = mean(correction, na.rm = TRUE))
count_novel_BRCA <- count_sel_gene %>%
  filter(grepl("BRCA",lib), grepl("unique", common_type)) %>%
  group_by(CL_lib) %>%
  summarise(n_novel = sum(count_pairs))
count_BRCA <- full_join(count_amp_BRCA, count_novel_BRCA) %>%
  replace(is.na(.), 0) %>%
  mutate(tissue = "BRCA")

count_tot <- rbind(count_COLO, count_BRCA)
pl <- ggplot(count_tot, aes(x = perc_amp, y = n_novel)) +
  geom_point(aes(color = tissue), size = 2) + 
  stat_smooth(method = "lm",
              formula = y ~ x, 
              se = FALSE) + 
  scale_color_manual(values = c("lightblue", "red")) + 
  theme_bw() +
  ylab("N. unique selected gene pairs") + 
  xlab("% of pairs amplified (Gene1 or Gene2)") 
pl
ggsave(pl,  
       filename = sprintf("%s/ALL_CLs/ALL_CLs_nUniqueGenePairs_VS_PercAmp.png", fold_output), width = 5, height = 4)

pl <- ggplot(count_tot, aes(x = mean_correction_dual, y = n_novel)) +
  geom_point(aes(color = tissue), size = 2) + 
  stat_smooth(method = "lm",
              formula = y ~ x, 
              se = FALSE) + 
  scale_color_manual(values = c("lightblue", "red")) + 
  theme_bw() +
  ylab("N. unique selected gene pairs") + 
  xlab("Mean correction") 
pl
ggsave(pl,  
       filename = sprintf("%s/ALL_CLs/ALL_CLs_nUniqueGenePairs_VS_meanCorrection.png", fold_output), width = 5, height = 4)

#### correction across libraries ####
# dual_FC_COLO <- output$dual_FC %>%
#   filter(lib == "COLO")
