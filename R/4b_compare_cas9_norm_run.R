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
library(VennDiagram)
library(argparse)
pdf(NULL) # avoid Rplots.pdf from being generated

parser <- ArgumentParser(description = "Compare the cas9Neg normalization")
parser$add_argument("--fold_input", type = "character", nargs = "*", help = "path to folder output, different cas9neg norm")
parser$add_argument("--name_input", type = "character", nargs = "*", help = "name cas9neg version")
parser$add_argument("--fold_output", type = "character", help = "path to folder to save output")

args <- parser$parse_args()
fold_input <- args$fold_input
name_input <- args$name_input
fold_output <- args$fold_output

########
# root_path <- "/group/iorio/lucia/"
# fold_output <- sprintf("%sCRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/BATCH_CORRECTED/SUMMARY_cas9neg_norm/", root_path)
# common <- sprintf("%sCRISPR_combinatorial/CRISPRcleanRatSquared/DATA_FREEZE_v4_NOV_2023/BATCH_CORRECTED/", root_path)
# fold_input <- c(sprintf("%scavg/ALL_CLs/", common), sprintf("%sc91/ALL_CLs/", common), sprintf("%sc92/ALL_CLs/", common))
# name_input <- c("cavg", "c91", "c92")
########

# plot model performances
print("Differences in model performances")
file <- lapply(fold_input, function(x) sprintf("%sALL_CLs_avg_model_perf_DoubleGenePairs.tsv", x))

model_perf <- mapply(function(x, y) readr::read_tsv(x) %>% 
                       dplyr::mutate(type = y), 
                     x = file, 
                     y = name_input, 
                     SIMPLIFY = FALSE) %>%
  bind_rows() %>%
  dplyr::mutate(CL_lib = paste(CL, lib, sep = "_"))

# plot
my_comparisons <- list( c("cavg", "c91"), c("cavg", "c92"))
pl <- ggboxplot(model_perf, x = "type", y = "mean_corr",
          color = "type", palette = "jco") + 
  stat_compare_means(comparisons = my_comparisons, 
                     paired = TRUE, 
                     method = "t.test") +
  geom_point(aes(color = type)) +
  xlab("") + 
  ylab("Mean model performance") +
  theme_bw() +
  theme(legend.position = "none")
pl
ggsave(filename = sprintf("%scompare_performances.pdf", fold_output), 
       plot = pl, width = 4, height = 4)

pl <- ggboxplot(model_perf, x = "type", y = "mean_correction_dual",
                color = "type", palette = "jco") + 
  stat_compare_means(comparisons = my_comparisons, 
                     paired = TRUE, 
                     method = "t.test") +
  geom_point(aes(color = type)) +
  xlab("") + 
  ylab("Mean correction dual") +
  theme_bw() +
  theme(legend.position = "none")
pl
ggsave(filename = sprintf("%scompare_correction_dual.pdf", fold_output), 
       plot = pl, width = 4, height = 4)

# plot gene selection
print("Differences in gene selection")
file <- lapply(fold_input, function(x) sprintf("%sselection_gene_letANDsyn.txt", x))

sel_output <- mapply(function(x, y) readr::read_tsv(x) %>% 
                       dplyr::mutate(cas9neg = y), 
                     x = file, 
                     y = name_input, 
                     SIMPLIFY = FALSE) %>%
  bind_rows() %>%
  dplyr::mutate(type_new = case_when(common == "in common" ~ "in common", 
                                     common == "unique" & type == "Uncorrected" ~ "unique (Uncorrected)", 
                                     common == "unique" & type == "Corrected" ~ "unique (CCR^2 Corrected)")) %>%
  dplyr::mutate(lib = factor(lib, levels = c("COLO", "BRCA")))

pl <- ggplot(sel_output, 
             aes(x = type_new, fill = cas9neg)) +
  geom_bar(position = position_dodge()) + 
  facet_wrap(.~lib, ncol = 1) + 
  theme_bw() + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  coord_flip() + 
  xlab("") + 
  ylab("N. selected gene pairs across all CLs") 
pl
ggsave(filename = sprintf("%scount_selected_gene_pairs.pdf", fold_output), 
       plot = pl, width = 5.2, height = 3.5)

pl <- ggplot(sel_output, 
             aes(x = CL, fill = cas9neg)) +
  geom_bar(position = position_dodge()) + 
  facet_wrap(.~lib, ncol = 1, scales = "free_y") + 
  theme_bw() + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  coord_flip() + 
  xlab("") + 
  ylab("N. selected gene pairs across all CLs") 
pl
ggsave(filename = sprintf("%scount_selected_gene_pairs_perCLs.pdf", fold_output), 
       plot = pl, width = 5.5, height = 7)

# intersection (all hits)
list_gene_pairs <- sel_output %>% 
  dplyr::mutate(cas9neg = factor(cas9neg, levels = name_input), 
                CL_lib_Gene_Pair = paste(CL_lib, Gene_Pair, sep = "_")) %>%
  dplyr::group_split(cas9neg) %>%
  lapply(function(x) x$CL_lib_Gene_Pair)

VennDiagram::venn.diagram(
  x = list_gene_pairs,
  category.names = name_input,
  filename = sprintf("%svenndiagram_selected_gene_pairs.png", fold_output),
  height = 620, 
  width = 620, 
  resolution = 300,
  output = TRUE
)



