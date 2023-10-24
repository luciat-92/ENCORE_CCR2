#### auxiliary functions ####
# combine top results
get_CNA <- function(copy_number_file, CL_name){
  
  # add CN info
  CNA <- readr::read_table(copy_number_file, 
                           show_col_types = FALSE) %>%
    dplyr::mutate(
      CHROM = dplyr::case_when(
        CHROM == "chrX" ~ "chr23",
        CHROM == "chrY" ~ "chr24",
        TRUE ~ as.character(CHROM))
      # .default = as.character(CHROM))
    )
  
  if ("Sampleid" %in% colnames(CNA)) {
    CL_name_unif <- toupper(str_replace_all(CL_name, "[-|_]", ""))
    CNA <- CNA %>% 
      dplyr::filter(toupper(str_replace_all(Sampleid, "[-|_]", "")) %in% CL_name_unif)
    if (nrow(CNA) == 0) {
      stop(sprintf("no CN info available for %s", CL_name))
    }
  }
  
  return(CNA)
  
}


add_CNA <- function(
    CNA, 
    dual_FC
) { 
  
  # use GenomicRanges package to match
  CNA <- CNA %>%
    dplyr::mutate(dplyr::across("CHROM", .fns = \(x) 
                                str_replace(x, pattern = "chr", replacement = "")))
  
  CNA_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    CNA[, c("CHROM", "start", "end")],
    seqnames.field = "CHROM",
    start.field = "start",
    end.field = "end")
  
  levels_cat <- c("High_loss", "Loss", "Neutral", "Soft_gain", "Amp")
  
  # remove pairs with NA chr (nontarget)
  dual_FC_no_nt <- dual_FC %>%
    dplyr::filter(!is.na(sgRNA1_Chr) & !is.na(sgRNA2_Chr))
  
  pos1_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    dual_FC_no_nt[, c("sgRNA1_Chr", "sgRNA1_Start", "sgRNA1_End")],
    seqnames.field = "sgRNA1_Chr",
    start.field = "sgRNA1_Start",
    end.field = "sgRNA1_End")
  
  pos2_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    dual_FC_no_nt[, c("sgRNA2_Chr", "sgRNA2_Start", "sgRNA2_End")],
    seqnames.field = "sgRNA2_Chr",
    start.field = "sgRNA2_Start",
    end.field = "sgRNA2_End")
  
  pos1_overlap <- GenomicRanges::findOverlaps(
    pos1_genom_range, 
    CNA_genom_range)
  
  pos2_overlap <- GenomicRanges::findOverlaps(
    pos2_genom_range, 
    CNA_genom_range)
  
  dual_FC_no_nt_withCNA <- dual_FC_no_nt %>%
    dplyr::mutate(Gene1_CN = NA, Gene1_CN_cat = NA,
                  Gene2_CN = NA, Gene2_CN_cat = NA) 
  dual_FC_no_nt_withCNA$Gene1_CN[queryHits(pos1_overlap)] <- CNA$C[subjectHits(pos1_overlap)]
  dual_FC_no_nt_withCNA$Gene2_CN[queryHits(pos2_overlap)] <- CNA$C[subjectHits(pos2_overlap)]
  dual_FC_no_nt_withCNA$Gene1_CN_cat[queryHits(pos1_overlap)] <- CNA$CN_category[subjectHits(pos1_overlap)]
  dual_FC_no_nt_withCNA$Gene2_CN_cat[queryHits(pos2_overlap)] <- CNA$CN_category[subjectHits(pos2_overlap)]
  
  dual_FC_no_nt_withCNA  <- dual_FC_no_nt_withCNA %>% 
    dplyr::mutate(
      Gene1_CN = round(Gene1_CN),
      Gene2_CN = round(Gene2_CN),
      Prod_CN = Gene1_CN * Gene2_CN, 
      Sum_CN = Gene1_CN + Gene2_CN, 
      Gene1_CN_cat = factor(Gene1_CN_cat, levels = levels_cat), 
      Gene2_CN_cat = factor(Gene2_CN_cat, levels = levels_cat), 
      Sum_CN_cat = as.numeric(Gene1_CN_cat) - 1 + as.numeric(Gene2_CN_cat) - 1) %>%
    dplyr::filter(!is.na(Prod_CN))
  
  dual_FC_no_nt_withCNA$Max_CN <- mapply(function(x, y) max(x,y), 
                                         x = dual_FC_no_nt_withCNA$Gene1_CN, 
                                         y = dual_FC_no_nt_withCNA$Gene2_CN, SIMPLIFY = T)
  dual_FC_no_nt_withCNA$Max_CN_cat <- mapply(function(x, y) max(x,y), 
                                             x = as.numeric(dual_FC_no_nt_withCNA$Gene1_CN_cat) - 1, 
                                             y = as.numeric(dual_FC_no_nt_withCNA$Gene2_CN_cat) - 1, SIMPLIFY = T)
  
  # consider nontarget-ess pairs
  dual_FC_nt <- dual_FC %>%
    dplyr::filter((is.na(sgRNA1_Chr) | is.na(sgRNA2_Chr)) & 
                    !(is.na(sgRNA1_Chr) & is.na(sgRNA2_Chr)))
  
  tmp1 <- dual_FC_nt %>% filter(!(is.na(sgRNA1_Chr)))
  pos1_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    tmp1[, c("sgRNA1_Chr", "sgRNA1_Start", "sgRNA1_End")],
    seqnames.field = "sgRNA1_Chr",
    start.field = "sgRNA1_Start",
    end.field = "sgRNA1_End")
  
  tmp2 <- dual_FC_nt %>% filter(!(is.na(sgRNA2_Chr)))
  pos2_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    tmp2[, c("sgRNA2_Chr", "sgRNA2_Start", "sgRNA2_End")],
    seqnames.field = "sgRNA2_Chr",
    start.field = "sgRNA2_Start",
    end.field = "sgRNA2_End")
  
  pos1_overlap <- GenomicRanges::findOverlaps(
    pos1_genom_range, 
    CNA_genom_range)
  
  pos2_overlap <- GenomicRanges::findOverlaps(
    pos2_genom_range, 
    CNA_genom_range)
  
  tmp1_withCNA <- tmp1 %>%
    dplyr::mutate(Gene1_CN = NA, Gene1_CN_cat = NA,
                  Gene2_CN = NA, Gene2_CN_cat = NA) 
  tmp1_withCNA$Gene1_CN[queryHits(pos1_overlap)] <- CNA$C[subjectHits(pos1_overlap)]
  tmp1_withCNA$Gene1_CN_cat[queryHits(pos1_overlap)] <- CNA$CN_category[subjectHits(pos1_overlap)]
  tmp1_withCNA <- tmp1_withCNA %>% 
    dplyr::mutate(
      Gene1_CN = round(Gene1_CN),
      Prod_CN = NA,
      Sum_CN = round(Gene1_CN),
      Gene1_CN_cat = factor(Gene1_CN_cat, levels = levels_cat), 
      Max_CN = round(Gene1_CN), 
      Max_CN_cat = as.numeric(Gene1_CN_cat) - 1) 
  
  tmp2_withCNA <- tmp2 %>%
    dplyr::mutate(Gene1_CN = NA, Gene1_CN_cat = NA,
                  Gene2_CN = NA, Gene2_CN_cat = NA) 
  tmp2_withCNA$Gene2_CN[queryHits(pos2_overlap)] <- CNA$C[subjectHits(pos2_overlap)]
  tmp2_withCNA$Gene2_CN_cat[queryHits(pos2_overlap)] <- CNA$CN_category[subjectHits(pos2_overlap)]
  tmp2_withCNA <- tmp2_withCNA %>% 
    dplyr::mutate(
      Gene2_CN = round(Gene2_CN),
      Prod_CN = NA,
      Sum_CN = round(Gene2_CN),
      Gene2_CN_cat = factor(Gene2_CN_cat, levels = levels_cat), 
      Max_CN = round(Gene2_CN), 
      Max_CN_cat = as.numeric(Gene2_CN_cat) - 1) 
  
  dual_FC_nt_withCNA <- bind_rows(tmp1_withCNA, tmp2_withCNA)
  
  # combine
  dual_FC_withCNA <- bind_rows(dual_FC_no_nt_withCNA, dual_FC_nt_withCNA)
  
  dual_FC_withCNA$Gene1_CN <- factor(dual_FC_withCNA$Gene1_CN, 
                                     levels = sort(unique(dual_FC_withCNA$Gene1_CN)))
  dual_FC_withCNA$Gene2_CN <- factor(dual_FC_withCNA$Gene2_CN, 
                                     levels = sort(unique(dual_FC_withCNA$Gene2_CN)))
  dual_FC_withCNA$Sum_CN <- factor(dual_FC_withCNA$Sum_CN, 
                                   levels = sort(unique(dual_FC_withCNA$Sum_CN)))
  dual_FC_withCNA$Prod_CN <- factor(dual_FC_withCNA$Prod_CN, 
                                    levels = sort(unique(dual_FC_withCNA$Prod_CN)))
  dual_FC_withCNA$Max_CN  <- factor(dual_FC_withCNA$Max_CN, 
                                    levels = sort(unique(dual_FC_withCNA$Max_CN)))  
  
  dual_FC_withCNA <- dual_FC_withCNA %>%
    filter(!is.na(Max_CN)) # remove those with NA in a position
  
  return(dual_FC_withCNA)
  
}

add_CNA_single <- function(
    CNA, 
    single_FC
){
  
  # match with CNA info
  CNA <- CNA %>%
    dplyr::mutate(dplyr::across("CHROM", .fns = \(x) 
                                str_replace(x, pattern = "chr", replacement = "")))
  
  CNA_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    CNA[, c("CHROM", "start", "end")],
    seqnames.field = "CHROM",
    start.field = "start",
    end.field = "end")
  
  single_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    single_FC[, c("CHR", "startp", "endp")],
    seqnames.field = "CHR",
    start.field = "startp",
    end.field = "endp")
  
  overlap <- GenomicRanges::findOverlaps(
    single_genom_range, 
    CNA_genom_range)
  
  levels_cat <- c("High_loss", "Loss", "Neutral", "Soft_gain", "Amp")
  
  FC_withCNA <- single_FC %>%
    dplyr::mutate(Gene_CN = NA, 
                  Gene_CN_cat = NA) 
  FC_withCNA$Gene_CN[queryHits(overlap)] <- CNA$C[subjectHits(overlap)]
  FC_withCNA$Gene_CN_cat[queryHits(overlap)] <- CNA$CN_category[subjectHits(overlap)]
  FC_withCNA <- FC_withCNA %>%
    dplyr::mutate(Gene_CN = round(Gene_CN), 
                  Gene_CN_cat = factor(Gene_CN_cat, levels = levels_cat))
  
  return(FC_withCNA)
  
}


combine_top <- function(top_corrected, top_uncorrected) {
  
  if (nrow(top_uncorrected) == 0 & nrow(top_corrected) > 0) {
    
    df <- data.frame(bliss_zscore = c(top_corrected$bliss_zscore_corrected),
                     order_id = 1:nrow(top_corrected),
                     ID = top_corrected$ID,
                     Gene_Pair = top_corrected$Gene_Pair,
                     logFC = top_corrected$correctedFC_scaled,
                     type = rep("Corrected", nrow(top_corrected)))
    df$common <- "unique"
  }else{
    if (nrow(top_uncorrected) > 0 & nrow(top_corrected) == 0) {
      df <- data.frame(bliss_zscore = c(top_uncorrected$bliss_zscore),
                       order_id = 1:nrow(top_uncorrected),
                       ID = top_uncorrected$ID,
                       Gene_Pair = top_uncorrected$Gene_Pair,
                       logFC = top_uncorrected$avgFC_scaled,
                       type = rep("Uncorrected", nrow(top_uncorrected)))
      df$common <- "unique"
    }else{
      
      if(nrow(top_uncorrected) == 0 & nrow(top_corrected) == 0){
        df <- NULL
      }else{
        
        common_id <- intersect(top_uncorrected$ID, top_corrected$ID)
        df <- data.frame(bliss_zscore = c(top_uncorrected$bliss_zscore,
                                          top_corrected$bliss_zscore_corrected),
                         order_id = c(1:nrow(top_uncorrected),
                                      1:nrow(top_corrected)),
                         ID = c(top_uncorrected$ID,
                                top_corrected$ID),
                         Gene_Pair = c(top_uncorrected$Gene_Pair,
                                       top_corrected$Gene_Pair),
                         logFC = c(top_uncorrected$avgFC_scaled,
                                   top_corrected$correctedFC_scaled),
                         type = c(rep("Uncorrected",nrow(top_uncorrected)),
                                  rep("Corrected", nrow(top_corrected))))
        
        df$common <- "unique"
        df$common[df$ID %in% common_id] <- "in common"
      }
    }
  }
  return(df)  
}

load_output <- function(fold, CL_name, lib_name, copy_number_file) {
  
  ccr2_output <- NULL
  file_name <- sprintf("%s%s/%s/%s_dual_correctedFCs.RData", 
                       fold, CL_name, lib_name, CL_name)
  if (file.exists(file_name)) {
    ccr2_output <- get(load(file_name))
    CNA <- get_CNA(copy_number_file = copy_number_file, 
                   CL_name = CL_name)
    
    ccr2_output$dual_with_CNA <- add_CNA(
      dual_FC = ccr2_output$dual, 
      CNA = CNA)
    
    ccr2_output$single_FC_withCNA <- add_CNA_single(
      single_FC = ccr2_output$single, 
      CNA = CNA)
    
    ccr2_output$selection <- combine_top(
      top_corrected = ccr2_output$top_corrected, 
      top_uncorrected = ccr2_output$top_uncorrected)
    
    ccr2_output$selection_gene <- combine_top(
      top_corrected = ccr2_output$top_gene_corrected, 
      top_uncorrected = ccr2_output$top_gene_uncorrected)
    
    ccr2_output$thr_correction <- data.frame(
      min_thr = ccr2_output$min_correction, 
      max_thr = ccr2_output$max_correction, 
      min_correction = min(ccr2_output$dual$correction, na.rm = T),
      max_correction = max(ccr2_output$dual$correction, na.rm = T)
    )
    
  }
  return(ccr2_output)
  
}

load_output_single <- function(fold, CL_name, lib_name) {
  
  single_output <- NULL
  file_name <- sprintf("%s%s/%s/%s_foldChanges.RData", 
                       fold, CL_name, lib_name, CL_name)
  
  if (file.exists(file_name)) {
    single_output <- get(load(file_name))
    single_output$avgFC_notcentered <- rowMeans(single_output[, -c(1:2), drop = F])
    single_output$avgFC <- single_output$avgFC_notcentered - median(single_output$avgFC_notcentered)
  }
  return(single_output)
  
}

get_all_CLs <- function(model_encore_table, fold, copy_number_file) {
  
  dual_FC <- dual_FC_withCNA <- single_FC_withCNA <- 
    selection_letANDsyn <-  selection_gene_letANDsyn <- 
    model_perf <- thr_correction <- list()
  
  for (idx in seq_len(nrow(model_encore_table))) {
    
    CL_name <- model_encore_table$model_name_CMP[idx]
    lib_name <- model_encore_table$lib[idx]
    print(paste(CL_name, lib_name, sep = "_"))
    
    out <- load_output(fold = fold, 
                       CL_name = CL_name,
                       lib_name = lib_name, 
                       copy_number_file = copy_number_file)
    
    if (!is.null(out)) {
      thr_correction[[idx]] <- out$thr_correction %>%
        mutate(
          CL = CL_name, 
          lib = lib_name,
          CL_lib = paste(CL_name, lib_name, sep = "_"))
        
      dual_FC[[idx]] <- out$dual %>% 
        mutate(
          CL = CL_name, 
          lib = lib_name,
          CL_lib = paste(CL_name, lib_name, sep = "_"))
      
      dual_FC_withCNA[[idx]] <- out$dual_with_CNA %>% 
        mutate(
          CL = CL_name, 
          lib = lib_name,
          CL_lib = paste(CL_name, lib_name, sep = "_"))
      
      single_FC_withCNA[[idx]] <- out$single_FC_withCNA %>% 
        mutate(
          CL = CL_name, 
          lib = lib_name,
          CL_lib = paste(CL_name, lib_name, sep = "_"))
      
      model_perf[[idx]] <- out$model_perf %>% 
        mutate(
          CL = CL_name, 
          lib = lib_name,
          CL_lib = paste(CL_name, lib_name, sep = "_"))
      
      if (!is.null(out$selection)) {
        selection_letANDsyn[[idx]] <- out$selection %>% 
          mutate(
            CL = CL_name, 
            lib = lib_name,
            CL_lib = paste(CL_name, lib_name, sep = "_"))
      }
      
      if (!is.null(out$selection_gene)) {
        selection_gene_letANDsyn[[idx]] <- out$selection_gene %>% 
          mutate(
            CL = CL_name, 
            lib = lib_name,
            CL_lib = paste(CL_name, lib_name, sep = "_"))
      }
      
      
    }else{
      print("No single screen available")
    }
  }
  
  levels_cat <- c("High_loss", "Loss", "Neutral", "Soft_gain", "Amp")
  
  dual_FC <- dplyr::bind_rows(dual_FC)
  dual_FC_withCNA <- dplyr::bind_rows(dual_FC_withCNA) %>%
    mutate(Gene1_CN = 
             factor(Gene1_CN, levels = sort(unique(as.numeric(as.character(Gene1_CN)))))) %>%
    mutate(Gene2_CN = 
             factor(Gene2_CN, levels = sort(unique(as.numeric(as.character(Gene2_CN)))))) %>%
    mutate(Sum_CN = 
             factor(Sum_CN, levels = sort(unique(as.numeric(as.character(Sum_CN)))))) %>%
    mutate(Prod_CN = 
             factor(Prod_CN, levels = sort(unique(as.numeric(as.character(Prod_CN)))))) %>%
    mutate(Max_CN = 
             factor(Max_CN, levels = sort(unique(as.numeric(as.character(Max_CN)))))) %>%
    mutate(Max_CN_cat = cut(Max_CN_cat, length(levels_cat), labels = levels_cat))
  
  single_FC_withCNA <- dplyr::bind_rows(single_FC_withCNA) %>%
    mutate(Gene_CN = 
             factor(Gene_CN, levels = sort(unique(as.numeric(as.character(Gene_CN)))))) 
  
  selection_letANDsyn <- dplyr::bind_rows(selection_letANDsyn)
  selection_gene_letANDsyn <- dplyr::bind_rows(selection_gene_letANDsyn)
  model_perf <- dplyr::bind_rows(model_perf)
  thr_correction <- dplyr::bind_rows(thr_correction)
  
  return(list(dual_FC = dual_FC, 
              dual_FC_withCNA = dual_FC_withCNA, 
              model_perf = model_perf,
              single_FC_withCNA = single_FC_withCNA, 
              selection_letANDsyn = selection_letANDsyn, 
              selection_gene_letANDsyn = selection_gene_letANDsyn, 
              thr_correction = thr_correction))
}

get_all_CLs_singleGW <- function(model_encore_table, fold) {
  
  single_FC <- list()
  CL_list <- model_encore_table[!duplicated(model_encore_table$model_name_CMP), ]
  for (idx in seq_len(nrow(CL_list))) {
    
    CL_name <- CL_list$model_name_CMP[idx]
    lib_name <- CL_list$lib[idx]
    print(CL_name)
    
    out <- load_output_single(fold = fold, 
                              CL_name = CL_name,
                              lib_name = lib_name)
    
    if (!is.null(out)) {
      single_FC[[idx]] <- out %>% 
        mutate(CL = CL_name)
    }else{
      print("No single screen available")
    }
  }
  
  single_FC <- dplyr::bind_rows(single_FC) 
  return(list(single_FC = single_FC))
}

# plots:
### mod to plot category
plotCNAdensity <- function(
    dual_FC_withCNA,
    saveToFig = F, 
    saveFormat = "pdf",
    outdir = "./", 
    EXPname = "", 
    CN_cat_thr, 
    essential_genes = NULL) {
  
  var_name <- "avgFC"
  var_name_corrected <- "correctedFC"
  xlab_name <- "logFC"
  
  # df <- dual_FC_withCNA %>%
  #   dplyr::filter(!is.na(correction)) %>%
  #   dplyr::mutate(CN_class = dplyr::case_when(
  #     as.numeric(as.character(Gene1_CN)) >= !!(CN_thr) | 
  #       as.numeric(as.character(Gene2_CN)) >= !!(CN_thr) ~ 
  #       sprintf("CN guide1 or guide2 >= %i", CN_thr), 
  #     TRUE ~ "Others")) 
  
  df <- dual_FC_withCNA %>%
    dplyr::filter(!is.na(correction)) %>%
    dplyr::mutate(CN_class = dplyr::case_when(
      Max_CN_cat >= CN_cat_thr ~ 
        sprintf("CN guide1 or guide2 >= %i", CN_cat_thr), 
      TRUE ~ "Others")) 
  
  df_plot <- data.frame(
    logFC = c(df[, var_name], df[, var_name_corrected]),
    CN_class = rep(df$CN_class, 2), 
    type = c(rep("pre-CRISPRCleanR^2", nrow(df)), 
             rep("post-CRISPRCleanR^2", nrow(df))), 
    Gene1 = rep(df$Gene1, 2), 
    Gene2 = rep(df$Gene2, 2)) 
  df_plot$type <- factor(df_plot$type, 
                         levels = c("pre-CRISPRCleanR^2", "post-CRISPRCleanR^2"))
  
  if (!is.null(essential_genes)) {
    df_plot <- df_plot %>% 
      filter(!(Gene1 %in% essential_genes | Gene2 %in% essential_genes))
    title_plot <- paste("Exclude Essential genes")
    add_to_name <- paste0("_rmEss")
  } else {
    title_plot = NULL
    add_to_name = ""
  }
  
  pl <- ggplot(df_plot, aes(x = logFC, fill = CN_class)) + 
    geom_density(alpha = 0.7) + 
    theme_bw() + 
    facet_wrap(.~type, nrow = 2) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 0.8) +
    theme(axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          legend.position = "bottom", legend.title = element_blank()) +
    scale_fill_manual(values = c("red", "grey70")) +
    xlab(xlab_name) +
    ggtitle(title_plot)
  
  print(pl)
  
  if (saveToFig) {
    file_name <- sprintf("%s%s_density%s_highCNcat%s.%s", 
                         outdir, EXPname, xlab_name, add_to_name, saveFormat)
    ggplot2::ggsave(filename =  file_name, plot = pl, width = 5, height = 5)
  }
  
}

plotCNA <- function(
    dual_FC_withCNA, 
    saveToFig = FALSE, 
    saveFormat = "pdf",
    outdir = "./", 
    EXPname ="", 
    essential_genes = NULL, 
    exclude_ess = TRUE) {
  
  
  var_name <- "avgFC"
  var_name_corrected <- "correctedFC"
  ylab_name <- "logFC"
  
  title_plot <- NULL
  add_to_name <- "" 
  height_pl <- 4.5
  if (!is.null(essential_genes) & exclude_ess) {
    dual_FC_withCNA <- dual_FC_withCNA %>% 
      filter(!(Gene1 %in% essential_genes | Gene2 %in% essential_genes))
    title_plot <- paste("Exclude Essential genes")
    add_to_name <- paste0("_rmEss")
    height_pl <- 5
  } else {
    if (!exclude_ess) {
      dual_FC_withCNA <- dual_FC_withCNA %>% 
        filter((Gene1 %in% essential_genes | Gene2 %in% essential_genes))
      title_plot <- paste("Essential genes")
      add_to_name <- paste0("_Ess")
      height_pl <- 5
    }
  }
  
  
  # remove outliers
  out_lim1 <- dual_FC_withCNA %>% 
    filter(!is.na(get(var_name))) %>%
    group_by(Max_CN_cat) %>% 
    summarise(I_low = max(min(get(var_name)), quantile(get(var_name), 0.25)  - 1.5*stats::IQR(get(var_name))), 
              I_up = min(max(get(var_name)), quantile(get(var_name), 0.75)  + 1.5*stats::IQR(get(var_name))))
  out_lim2 <- dual_FC_withCNA %>% 
    filter(!is.na(get(var_name_corrected))) %>%
    group_by(Max_CN_cat) %>% 
    summarise(I_low = max(min(get(var_name_corrected)), quantile(get(var_name_corrected), 0.25)  - 1.5*stats::IQR(get(var_name_corrected))), 
              I_up = min(max(get(var_name_corrected)), quantile(get(var_name_corrected), 0.75)  + 1.5*stats::IQR(get(var_name_corrected))))
  lim_low <- min(c(out_lim1$I_low, out_lim2$I_low))
  lim_up <- max(c(out_lim1$I_up, out_lim2$I_up))
  
  df_plot <- data.frame(
    CL_lib = rep(dual_FC_withCNA$CL_lib, 2),
    info_subtype = rep(dual_FC_withCNA$info_subtype, 2),
    logFC = c(dual_FC_withCNA[, var_name, drop = T], dual_FC_withCNA[, var_name_corrected, drop = T]), 
    Sum_CN = rep(dual_FC_withCNA$Sum_CN, 2), 
    Gene1_CN = rep(dual_FC_withCNA$Gene1_CN, 2), 
    Gene2_CN = rep(dual_FC_withCNA$Gene2_CN, 2), 
    Gene1_CN_cat = rep(dual_FC_withCNA$Gene1_CN_cat, 2), 
    Gene2_CN_cat = rep(dual_FC_withCNA$Gene2_CN_cat, 2), 
    Prod_CN = rep(dual_FC_withCNA$Prod_CN, 2), 
    Max_CN = rep(dual_FC_withCNA$Max_CN, 2), 
    Max_CN_cat = rep(dual_FC_withCNA$Max_CN_cat, 2), 
    type = c(rep("pre-CRISPRCleanR^2", nrow(dual_FC_withCNA)), 
             rep("post-CRISPRCleanR^2", nrow(dual_FC_withCNA))), 
    Gene1 = rep(dual_FC_withCNA$Gene1, 2), 
    Gene2 = rep(dual_FC_withCNA$Gene2, 2))
  df_plot$type <- factor(df_plot$type, 
                         levels = c("pre-CRISPRCleanR^2", "post-CRISPRCleanR^2"))
  df_plot$comb_CN <- sprintf("%s ~ %s", df_plot$Gene1_CN, df_plot$Gene2_CN)
  
  ### mod to plot only CAT ###
  pl_CN <- ggplot(df_plot, aes(x = Max_CN_cat, y = logFC, 
                               fill = type)) + 
    geom_boxplot(outlier.size = 0.1, outlier.colour = "transparent") + 
    ggpubr::stat_compare_means(aes(group = type), 
                               label = "p.signif", 
                               method = "wilcox.test",
                               paired = TRUE, 
                               label.y = lim_up - 0.05) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme_bw() + 
    theme(axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "bottom") +
    scale_fill_brewer(palette = "Paired") +
    xlab("Max CN guide1 & guide2") + 
    ylab(ylab_name) + 
    ggtitle(title_plot) +
    coord_cartesian(ylim = c(lim_low, lim_up))
  print(pl_CN)
  
  df_plot <- df_plot %>%
    filter(!(is.na(Gene1_CN) | is.na(Gene2_CN))) # remove any pairs with NA
  
  df_plot$Gene2_CN <- sprintf("CN guide2: %s", df_plot$Gene2_CN)
  df_plot$Gene2_CN <- factor(
    df_plot$Gene2_CN, 
    levels = sprintf("CN guide2: %s", sort(as.numeric(levels(dual_FC_withCNA$Gene2_CN)))))
  
  df_plot$Gene2_CN_cat <- sprintf("CN guide2: %s", df_plot$Gene2_CN_cat)
  df_plot$Gene2_CN_cat <- factor(
    df_plot$Gene2_CN_cat, 
    levels = sprintf("CN guide2: %s", levels(dual_FC_withCNA$Gene2_CN_cat)))
  
  out_lim1 <- dual_FC_withCNA %>% 
    filter(!is.na(get(var_name))) %>%
    group_by(Gene1_CN_cat, Gene2_CN_cat) %>% 
    summarise(I_low = max(min(get(var_name)), quantile(get(var_name), 0.25)  - 1.5*stats::IQR(get(var_name))), 
              I_up = min(max(get(var_name)), quantile(get(var_name), 0.75)  + 1.5*stats::IQR(get(var_name))))
  
  out_lim2 <- dual_FC_withCNA %>% 
    filter(!is.na(get(var_name_corrected))) %>%
    group_by(Gene1_CN_cat, Gene2_CN_cat) %>% 
    summarise(I_low = max(min(get(var_name_corrected)), quantile(get(var_name_corrected), 0.25)  - 1.5*stats::IQR(get(var_name_corrected))), 
              I_up = min(max(get(var_name_corrected)), quantile(get(var_name_corrected), 0.75)  + 1.5*stats::IQR(get(var_name_corrected))))
  lim_low <- min(c(out_lim1$I_low, out_lim2$I_low))
  lim_up <- max(c(out_lim1$I_up, out_lim2$I_up))
  
  pl_CN_comb <- ggplot(df_plot, aes(x = Gene1_CN_cat, y = logFC, fill = type)) + 
    geom_boxplot(outlier.size = 0.1, outlier.colour = "transparent") + 
    theme_bw() + 
    facet_wrap(.~Gene2_CN_cat, nrow = 1) +
    ggpubr::stat_compare_means(aes(group = type), 
                               label = "p.signif", 
                               method = "wilcox.test",
                               paired = TRUE, 
                               label.y = lim_up - 0.05) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          legend.position = "bottom", 
          legend.title = element_blank(), 
          axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_brewer(palette = "Paired") +
    xlab("CN guide1") + 
    ylab(ylab_name) + 
    ggtitle(title_plot) + 
    coord_cartesian(ylim = c(lim_low, lim_up))
  print(pl_CN_comb)
  
  if (saveToFig) {
    file_name_maxCN <- sprintf("%s%s_MaxCNcat_vs_%s%s.%s", 
                               outdir, EXPname, ylab_name, add_to_name, saveFormat)
    file_name_combCN <- sprintf("%s%s_CombCNcat_vs_%s%s.%s", 
                                outdir, EXPname, ylab_name, add_to_name, saveFormat)
    ggsave(filename = file_name_maxCN, plot = pl_CN, width = 6, height = height_pl)
    ggsave(filename = file_name_combCN, plot = pl_CN_comb, width = 14, height = height_pl)
  }
  
}

plotCNA_diff <- function(
    dual_FC_withCNA, 
    saveToFig = FALSE, 
    saveFormat = "pdf",
    outdir = "./", 
    EXPname ="", 
    essential_genes = NULL, 
    exclude_ess = TRUE) {
  
  ylab_name <- "mean(logFC post - logFC pre)"
  df_plot <- dual_FC_withCNA
  
  title_plot <- NULL
  add_to_name <- "" 
  height_pl <- 4
  if (!is.null(essential_genes) & exclude_ess) {
    df_plot <- df_plot %>% 
      filter(!(Gene1 %in% essential_genes | Gene2 %in% essential_genes))
    title_plot <- paste("Exclude Essential genes")
    add_to_name <- paste0("_rmEss")
    height_pl <- 4.5
  } else {
    if (!exclude_ess) {
      df_plot <- df_plot %>% 
        filter((Gene1 %in% essential_genes | Gene2 %in% essential_genes))
      title_plot <- paste("Essential genes")
      add_to_name <- paste0("_Ess")
      height_pl <- 4.5
    }
  }
  
  ### mod to plot only CAT ###
  pl_CN <- ggplot(df_plot, aes(x = Max_CN_cat, y = correction)) + 
    # geom_violin() + 
    stat_summary(fun = mean,
                 geom = "pointrange",
                 fun.min = function(x) mean(x) - sd(x),
                 fun.max = function(x) mean(x) + sd(x)) + 
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme_bw() + 
    theme(axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "bottom") +
    xlab("Max CN guide1 & guide2") + 
    ylab(ylab_name) + 
    ggtitle(title_plot)
  print(pl_CN)
  
  df_plot <- df_plot %>%
    dplyr::filter(!(is.na(Gene1_CN) | is.na(Gene2_CN)))
  
  df_plot$Gene2_CN <- sprintf("CN guide2: %s", df_plot$Gene2_CN)
  df_plot$Gene2_CN <- factor(
    df_plot$Gene2_CN, 
    levels = sprintf("CN guide2: %s", sort(as.numeric(levels(dual_FC_withCNA$Gene2_CN)))))
  
  df_plot$Gene2_CN_cat <- sprintf("CN guide2: %s", df_plot$Gene2_CN_cat)
  df_plot$Gene2_CN_cat <- factor(
    df_plot$Gene2_CN_cat, 
    levels = sprintf("CN guide2: %s", levels(dual_FC_withCNA$Gene2_CN_cat)))
  
  pl_CN_comb <- ggplot(df_plot, aes(x = Gene1_CN_cat, y = correction)) + 
    stat_summary(fun = mean,
                 geom = "pointrange",
                 fun.min = function(x) mean(x) - sd(x),
                 fun.max = function(x) mean(x) + sd(x)) + 
    theme_bw() + 
    facet_wrap(.~Gene2_CN_cat, nrow = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          legend.position = "bottom", 
          legend.title = element_blank(), 
          axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_brewer(palette = "Paired") +
    xlab("CN guide1") + 
    ylab(ylab_name) + 
    ggtitle(title_plot)
  
  print(pl_CN_comb)
  
  if (saveToFig) {
    file_name_maxCN <- sprintf("%s%s_MaxCNcat_vs_avgCorr%s.%s", 
                               outdir, EXPname,  add_to_name, saveFormat)
    file_name_combCN <- sprintf("%s%s_CombCNcat_vs_avgCorr%s.%s", 
                                outdir, EXPname,  add_to_name, saveFormat)
    ggsave(filename = file_name_maxCN, plot = pl_CN, width = 6, height = height_pl)
    ggsave(filename = file_name_combCN, plot = pl_CN_comb, width = 14, height = height_pl)
  }
  
}

plotCNA_count <- function(
    dual_FC_withCNA, 
    saveToFig = FALSE,
    essential_genes = NULL,
    saveFormat = "pdf",
    outdir = "./", 
    EXPname = "") {
  
  
  dual_FC_withCNA <-  dual_FC_withCNA %>%
    dplyr::mutate(type = dplyr::case_when(
      (Gene1 %in% essential_genes | Gene2 %in% essential_genes) ~ "Gene1 or Gene2 Core Essential",
      TRUE ~ "Otherwise"))
  
  df_plot_max <- dual_FC_withCNA %>%
    dplyr::group_by(Max_CN_cat, type) %>%
    dplyr::summarise(count = n(), 
                     CL = length(unique(CL)), 
                     Gene1_n = length(unique(Gene1)), 
                     Gene2_n = length(unique(Gene2))) %>%
    dplyr::ungroup() 
  
  df_plot_comb <- dual_FC_withCNA %>%
    dplyr::filter(!(is.na(Gene1_CN) | is.na(Gene2_CN))) %>%
    dplyr::group_by(Gene2_CN_cat, Gene1_CN_cat, type) %>%
    dplyr::summarise(count = n(), 
                     CL = length(unique(CL)), 
                     Gene1_n = length(unique(Gene1)), 
                     Gene2_n = length(unique(Gene2))) %>%
    dplyr::ungroup()
  df_plot_comb$Gene2_CN_cat <- sprintf("CN guide2: %s", df_plot_comb$Gene2_CN_cat)
  df_plot_comb$Gene2_CN_cat <- factor(
    df_plot_comb$Gene2_CN_cat, 
    levels = sprintf("CN guide2: %s", levels(dual_FC_withCNA$Gene2_CN_cat)))
  
  pl_count_max <- ggplot(df_plot_max, aes(x = Max_CN_cat, y = count, 
                                          fill = type)) + 
    geom_bar(stat = "identity", width = 0.5, color = "black", position = "dodge") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "top") +
    scale_fill_manual(values = c("red", "grey80")) +
    xlab("Max CN guide1 & guide2") + 
    ylab("N. of pairs") + 
    scale_y_continuous(trans = "log10", labels = scales::scientific_format())
  print(pl_count_max)
  
  pl_count_comb <- ggplot(df_plot_comb, aes(x = Gene1_CN_cat, y = count, fill = type)) + 
    geom_bar(stat = "identity", width = 0.5, color = "black", position = "dodge") + 
    theme_bw() + 
    facet_wrap(.~Gene2_CN_cat, nrow = 1) +
    theme(axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          legend.position = "top", 
          legend.title = element_blank(), 
          axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_manual(values = c("red", "grey80")) +
    xlab("CN guide1") + 
    ylab("N. of pairs") +
    scale_y_continuous(trans = "log10", labels = scales::scientific_format())
  print(pl_count_comb)
  
  
  if (saveToFig) {
    file_name_maxCN <- sprintf("%s%s_MaxCNcat_vs_%s.%s", 
                               outdir, EXPname, "count", saveFormat)
    file_name_combCN <- sprintf("%s%s_CombCNcat_vs_%s.%s", 
                                outdir, EXPname, "count", saveFormat)
    ggsave(filename = file_name_maxCN, plot = pl_count_max, width = 6, height = 4.5)
    ggsave(filename = file_name_combCN, plot = pl_count_comb, width = 14, height = 4.5)
  }
}

plot_model_perf <- function(model_perf, 
                            saveToFig = FALSE, 
                            saveFormat = "pdf",
                            outdir = "./", 
                            EXPname = ""){
  
  prepare_data <- function(tissue, model_perf){
    
    perf <- model_perf %>%
      filter(grepl(tissue, lib)) %>%
      mutate(type_position = sprintf("%s (%s)", type, position)) %>%
      mutate(type_position = factor(type_position , 
                                    levels = c("NONTARGET_PAIR (guide2)", "NONTARGET_PAIR (guide1)", "DOUBLE_CUT_PAIR (guide2)", "DOUBLE_CUT_PAIR (guide1)")))
    
    CL_ord <- perf %>% 
      filter(type == "DOUBLE_CUT_PAIR") %>%
      group_by(CL, lib) %>%
      summarise(sum_cor = sum(cor)) %>%
      ungroup() %>%
      group_by(CL) %>%
      summarise(avg_sum_cor = mean(sum_cor)) %>%
      ungroup() %>%
      arrange(avg_sum_cor) %>%
      pull(CL)
    
    pl1 <- ggplot(subset(perf, grepl("DOUBLE_CUT_PAIR",as.character(type_position))), 
                  aes(x = factor(CL, levels = CL_ord), y = cor, fill = position)) + 
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      facet_wrap(.~lib, nrow = 1) +
      theme_bw() + 
      theme(legend.title = element_blank(),
            legend.position = "bottom", 
            axis.text = element_text(size = 11), 
            axis.title = element_text(size = 12), 
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab("Pears. Corr.\nSingle VS Pseudo single") + 
      xlab("") +
      ggtitle("Double Gene Pairs") +
      coord_flip()
    
    pl2 <- ggplot(subset(perf, grepl("NONTARGET_PAIR", as.character(type_position))), 
                  aes(x = factor(CL, levels = CL_ord), y = cor, fill = position)) + 
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      facet_wrap(.~lib, nrow = 1) +
      theme_bw() + 
      theme(legend.title = element_blank(),
            legend.position = "bottom", 
            axis.text = element_text(size = 11), 
            axis.title = element_text(size = 12), 
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab("Pears. Corr.\nSingle VS Pseudo single") + 
      xlab("") +
      ggtitle("NON-TARGET Pairs") +
      coord_flip()
    
    pl <- ggarrange(plotlist = list(pl1, pl2), nrow = 1)
    return(pl)
    
  }
  
  # plot model performances
  pl1 <- prepare_data(tissue = "COLO", model_perf = model_perf)
  pl2 <- prepare_data(tissue = "BRCA", model_perf = model_perf)
  
  pl <- ggarrange(plotlist = list(pl1, pl2), 
                  ncol = 1, 
                  legend = "right",
                  common.legend = TRUE, 
                  align = "hv", 
                  heights = c(1, 0.7))
  print(pl)
  
  if (saveToFig) {
    file_name <- sprintf("%s%s_model_perf.%s", 
                         outdir, EXPname, saveFormat)
    ggsave(filename = file_name, plot = pl, width = 9, height = 8)
  }
}


plot_model_perf_vs_corr <- function(
    tissue_name,
    model_perf, 
    dual_FC_withCNA, 
    single_FC_withCNA,
    single_FC_gw,
    saveToFig = FALSE, 
    saveFormat = "pdf",
    outdir = "./", 
    EXPname = ""){
  
  prepare_data <- function(tissue, model_perf, dual_FC_withCNA, 
                           single_FC_withCNA, single_FC_gw, type){
    
    perf <- model_perf %>%
      filter(grepl(tissue, lib)) %>%
      filter(type == !!(type)) %>%
      group_by(CL, lib) %>%
      summarise(mean_corr = mean(cor))
    
    if (type == "DOUBLE_CUT_PAIR") {  
      count_amp <- dual_FC_withCNA %>%
        filter(!grepl("NONTARGET", info_subtype)) 
      
      dual <- dual_FC_withCNA %>%
        filter(!grepl("NONTARGET", info_subtype)) 
      
    }else{
      count_amp <- dual_FC_withCNA %>%
        filter(grepl("NONTARGET", info_subtype)) 
      
      dual <- dual_FC_withCNA %>%
        filter(grepl("NONTARGET", info_subtype)) 
    }
    
    single <- single_FC_withCNA %>%
      filter(genes %in% c(dual$Gene1, dual$Gene2))
    
    count_amp <- count_amp %>%
      filter(grepl(tissue, lib)) %>%
      group_by(CL, lib) %>%
      summarise(perc_pair_amp = sum(as.character(Max_CN_cat) %in% c("Amp"))/length(Max_CN_cat), 
                n_pair_amp = sum(as.character(Max_CN_cat) %in% c("Amp")))
    
    dual <- dual %>%
      filter(grepl(tissue, lib)) %>%
      group_by(CL, lib) %>%
      summarise(median_correction = median(correction, na.rm = TRUE), 
                mean_correction = mean(correction, na.rm = TRUE), 
                sd_correction = sd(correction, na.rm = TRUE))
    
    single <- single %>%
      filter(grepl(tissue, lib)) %>%
      group_by(CL, lib) %>%
      summarise(median_correction = median(correction, na.rm = TRUE), 
                mean_correction = mean(correction, na.rm = TRUE), 
                sd_correction = sd(correction, na.rm = TRUE))
    
    single_gw <- single_FC_gw %>% 
      filter(CL %in% dual$CL) %>%
      group_by(CL) %>%
      summarise(median_avgFC_singleGW = median(avgFC, na.rm = TRUE), 
                mean_avgFC_singleGW = mean(avgFC, na.rm = TRUE), 
                sd_avgFC_singleGW = sd(avgFC, na.rm = TRUE))
    
    tot <- full_join(dual, single, 
                     suffix = c("_dual", "_single"), 
                     by = c("CL", "lib"))
    tot <- full_join(tot, perf)
    tot <- full_join(tot, count_amp)
    tot <- full_join(tot, single_gw, by = c("CL"))
    
    return(tot)
  }
  
  # plot model performances
  df_gene <- prepare_data(tissue = tissue_name, 
                          model_perf = model_perf, 
                          dual_FC_withCNA = dual_FC_withCNA, 
                          single_FC_withCNA = single_FC_withCNA, 
                          single_FC_gw = single_FC_gw,
                          type = "DOUBLE_CUT_PAIR") %>%
    mutate(diff_mean = mean_correction_single - mean_correction_dual) %>%
    mutate(ratio_mean = mean_correction_dual/mean_correction_single) 
  
  df_nt <- prepare_data(tissue = tissue_name, 
                          model_perf = model_perf, 
                          dual_FC_withCNA = dual_FC_withCNA, 
                          single_FC_withCNA = single_FC_withCNA, 
                          single_FC_gw = single_FC_gw,
                          type = "NONTARGET_PAIR") %>%
    mutate(diff_mean = mean_correction_single - mean_correction_dual) %>%
    mutate(ratio_mean = mean_correction_dual/mean_correction_single) 
  
  plot_tmp <- function(df, df_nt, var_to_plot, ytitle, tissue_name){
    
    pl1 <- ggplot(df_nt, 
                  aes(x = mean_corr, 
                      y = get(var_to_plot))) + 
      geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") + 
      stat_cor(method = "pearson") + 
      geom_point(aes(color = perc_pair_amp), size = 2) +
      theme_bw() + 
      ylab(ytitle) + 
      xlab("Mean model performance") +
      theme(axis.text = element_text(size = 11), 
            axis.title = element_text(size = 12)) +
      scale_colour_gradient2(low = "red", mid = "grey50",high = "red", name = "% of pairs amplified") + 
      # scale_size(name = "% of pairs amplified") + 
      ggtitle(sprintf("NON-TARGET Pairs (%s)", tissue_name))
    
    pl2 <- ggplot(df, 
                  aes(x = mean_corr, 
                      y = get(var_to_plot))) + 
      geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") + 
      stat_cor(method = "pearson") + 
      geom_point(aes(color = perc_pair_amp), size = 2) +
      theme_bw() + 
      ylab(ytitle) + 
      xlab("Mean model performance") +
      theme(axis.text = element_text(size = 11), 
            axis.title = element_text(size = 12)) +
      scale_colour_gradient2(low = "red", mid = "grey50",high = "red", name = "% of pairs amplified") + 
      # scale_size(name = "% of pairs amplified") + 
      ggtitle(sprintf("Double Gene Pairs (%s)", tissue_name))
    
    pl <- ggarrange(plotlist = list(pl1, pl2), nrow = 1)
    
    return(pl)
  }
  
  pl1 <- plot_tmp(df_gene, 
                  df_nt,
                  var_to_plot = "mean_correction_dual", 
                  ytitle = "Mean correction dual", 
                  tissue_name = tissue_name)
  print(pl1)
  pl2 <- plot_tmp(
    df_gene, 
    df_nt,
    var_to_plot = "mean_correction_single", 
    ytitle = "Mean correction single (matching)", 
    tissue_name = tissue_name)
  print(pl2)
  
  if (saveToFig) {
    file_name <- sprintf("%s%s_%s_model_perf_VS_mean_correction_dual.%s", 
                         outdir, EXPname, tissue_name, saveFormat)
    ggsave(filename = file_name, plot = pl1, width = 10, height = 4)
    
    file_name <- sprintf("%s%s_%s_model_perf_VS_mean_correction_single.%s", 
                         outdir, EXPname, tissue_name, saveFormat)
    ggsave(filename = file_name, plot = pl2, width = 10, height = 4)
  }
  
  return(list(gene_pairs = df_gene, gene_nontarget = df_nt))
  
}

plot_count_amp <- function(dual_FC_withCNA, 
                           saveToFig = FALSE, 
                           saveFormat = "pdf",
                           outdir = "./", 
                           EXPname = ""){
  
  count_amp_COLO <- dual_FC_withCNA %>%
    filter(grepl("COLO", lib)) %>%
    group_by(CL, lib) %>%
    summarise(n_amp = sum(as.character(Max_CN_cat) == "Amp"))
  
  count_amp_BRCA <- dual_FC_withCNA %>%
    filter(grepl("BRCA", lib)) %>%
    group_by(CL, lib) %>%
    summarise(n_amp = sum(as.character(Max_CN_cat) == "Amp"))
  
  pl1 <- ggplot(count_amp_COLO, 
                aes(x = reorder(CL, n_amp), y = n_amp)) + 
    geom_bar(stat = "identity", width = 0.5) +
    facet_wrap(.~lib, nrow = 1) +
    theme_bw() + 
    ylab("N. Pairs with Gene1 or Gene2 Amp") + 
    xlab("") +
    theme(axis.text = element_text(size = 11), 
          axis.title = element_text(size = 12), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
  
  pl2 <- ggplot(count_amp_BRCA, 
                aes(x = reorder(CL, n_amp), y = n_amp)) + 
    geom_bar(stat = "identity", width = 0.5) +
    facet_wrap(.~lib, nrow = 1) +
    theme_bw() + 
    ylab("N. Pairs with Gene1 or Gene2 Amp") + 
    xlab("") +
    theme(axis.text = element_text(size = 11), 
          axis.title = element_text(size = 12), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
  
  pl <- ggarrange(plotlist = list(pl1, pl2), ncol = 1, 
                  common.legend = TRUE, 
                  align = "v", 
                  heights = c(1, 0.7))
  print(pl)
  
  if (saveToFig) {
    file_name <- sprintf("%s%s_count_pairs_amp.%s", 
                         outdir, EXPname, saveFormat)
    ggsave(filename = file_name, plot = pl, width = 5, height = 8)
  }
  
  return(count_amp = rbind(count_amp_COLO, count_amp_BRCA))
  
}

plot_CN_vs_logFC <- function(
    dual_FC_withCNA, 
    saveToFig = FALSE, 
    saveFormat = "pdf",
    outdir = "./", 
    EXPname ="", 
    essential_genes = NULL, 
    exclude_ess = TRUE) {
  
  var_name <- "avgFC"
  ylab_name <- "logFC"
  
  df_plot <- dual_FC_withCNA %>%
    rename(logFC = !!(var_name)) %>%
    filter(!is.na(logFC))
  
  title_plot <- NULL
  add_to_name <- "" 
  height_pl <- 5
  if (!is.null(essential_genes) & exclude_ess) {
    df_plot <- df_plot %>% 
      filter(!(Gene1 %in% essential_genes | Gene2 %in% essential_genes))
    title_plot <- paste("Exclude Essential genes")
    add_to_name <- paste0("_rmEss")
    height_pl <- 5.5
  } else {
    if (!exclude_ess) {
      df_plot <- df_plot %>% 
        filter((Gene1 %in% essential_genes | Gene2 %in% essential_genes))
      title_plot <- paste("Essential genes")
      add_to_name <- paste0("_Ess")
      height_pl <- 5.5
    }
  }
  
  ### mod to plot only CAT ###
  # remove outliers
  out_lim <- dual_FC_withCNA %>% 
    filter(!is.na(get(var_name))) %>%
    filter(!is.na(Max_CN_cat)) %>%
    group_by(Max_CN_cat) %>% 
    summarise(I_low = quantile(get(var_name), 0.25) - 1.5*(quantile(get(var_name), 0.75) - quantile(get(var_name), 0.25)), 
              I_up = quantile(get(var_name), 0.75) + 1.5*(quantile(get(var_name), 0.75) - quantile(get(var_name), 0.25)))
  lim_low <- min(out_lim$I_low)
  lim_up <- max(out_lim$I_up)
  
  pl_CN <- ggplot(df_plot, aes(x = Max_CN_cat, y = logFC)) + 
    geom_boxplot(outlier.size = 0.5) + 
    geom_hline(yintercept = 0, 
               color = "red", 
               linetype = "dashed", 
               size = 0.5) +
    theme_bw() + 
    theme(axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "bottom") +
    xlab("Max CN guide1 & guide2") + 
    ylab(ylab_name) + 
    ggtitle(title_plot) +
    coord_cartesian(ylim = c(lim_low, lim_up))
  print(pl_CN)
  
  df_plot <- df_plot %>%
    filter(!(is.na(Gene1_CN) | is.na(Gene2_CN))) # remove any pairs with NA
  
  df_plot$Gene2_CN <- sprintf("CN guide2: %s", df_plot$Gene2_CN)
  df_plot$Gene2_CN <- factor(
    df_plot$Gene2_CN, 
    levels = sprintf("CN guide2: %s", sort(as.numeric(levels(dual_FC_withCNA$Gene2_CN)))))
  
  df_plot$Gene2_CN_cat <- sprintf("CN guide2: %s", df_plot$Gene2_CN_cat)
  df_plot$Gene2_CN_cat <- factor(
    df_plot$Gene2_CN_cat, 
    levels = sprintf("CN guide2: %s", levels(dual_FC_withCNA$Gene2_CN_cat)))
  
  out_lim <- dual_FC_withCNA %>% 
    filter(!is.na(get(var_name))) %>%
    filter(!(is.na(Gene1_CN) | is.na(Gene2_CN))) %>%
    group_by(Gene1_CN_cat, Gene2_CN_cat) %>% 
    summarise(I_low = quantile(get(var_name), 0.25)  - 1.5*(quantile(get(var_name), 0.75) - quantile(get(var_name), 0.25)), 
              I_up = quantile(get(var_name), 0.75)  + 1.5*(quantile(get(var_name), 0.75) - quantile(get(var_name), 0.25))) %>%
    drop_na()
  lim_low <- min(out_lim$I_low)
  lim_up <- max(out_lim$I_up)
  
  pl_CN_comb <- ggplot(df_plot, aes(x = Gene1_CN_cat, y = logFC)) + 
    geom_boxplot(outlier.size = 0.3) + 
    theme_bw() + 
    facet_wrap(.~Gene2_CN_cat, nrow = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    theme(axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          legend.position = "bottom", 
          legend.title = element_blank(), 
          axis.text.x = element_text(angle = 60, hjust = 1)) +
    xlab("CN guide1") + 
    ylab(ylab_name) + 
    coord_cartesian(ylim = c(lim_low, lim_up)) +
    ggtitle(title_plot)
  print(pl_CN_comb)
  
  if (saveToFig) {
    file_name_maxCN <- sprintf("%s%s_pre_MaxCNcat_vs_%s%s.%s", 
                               outdir, EXPname, ylab_name, add_to_name, saveFormat)
    file_name_combCN <- sprintf("%s%s_pre_CombCNcat_vs_%s%s.%s", 
                                outdir, EXPname, ylab_name, add_to_name, saveFormat)
    ggsave(filename = file_name_maxCN, plot = pl_CN, width = 4, height = height_pl)
    ggsave(filename = file_name_combCN, plot = pl_CN_comb, width = 10, height = height_pl)
  }
}
