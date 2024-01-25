# compute average per cell line
compute_avglogFC_repetitions <- function(data_logFC){
  
  CL_names <- sort(unique(str_split_fixed(string = colnames(data_logFC), 
                                          pattern = "[_]", n = 2)[,1]))
  
  data_mean_logFC <- matrix(nrow = nrow(data_logFC), ncol = length(CL_names))
  
  for (idx in 1:length(CL_names)) {
    CL_id <- CL_names[idx]
    tmp <- data_logFC[,grepl(CL_id, colnames(data_logFC))]
    if (ncol(tmp) > 1) {
      tmp_mean <- rowMeans(tmp)
    }else{
      tmp_mean <- tmp[,1]
    }
    data_mean_logFC[,idx] <- tmp_mean
  }
  rownames(data_mean_logFC) <- rownames(data_logFC)
  colnames(data_mean_logFC) <- CL_names
  return(data_mean_logFC)
  
}

# load data (rds format)
load_data_rds <- function(fold_input, fold_lib, negcontrol){
  
  lib_name <- c("COLO1", "COLO2", "COLO3", 
                "BRCA1", "BRCA2", "BRCA3") %>% tolower()
  
  
  data <- library <- list()
  for (i in 1:length(lib_name)) {
    
    curr_lib <- lib_name[i]
    print(curr_lib)
    
    data[[i]] <- read_rds(sprintf("%s%s_%s.rds", fold_input, curr_lib, negcontrol))  
    # compute mean per same cell line across repetitions
    data[[i]] <- compute_avglogFC_repetitions(data[[i]])
    
    data_lib <- read_rds(sprintf("%sIDs_%s.rds", fold_input, curr_lib)) %>% as.data.frame()
    data_lib <- data_lib %>% arrange(ID)
    if (!identical(data_lib$ID, rownames(data[[i]]))) {
      print("reorder data_lib")
      data_lib <- data_lib[match(rownames(data[[i]]), data_lib$ID), ]
    }
    
    if (!"Gene" %in% colnames(data_lib)) {
      data_lib <- data_lib %>% 
        mutate(Gene_Pair = paste(Gene1, Gene2, sep = "~"), .before = "Gene1")
    }
    data[[i]] <- bind_cols(data_lib, as.data.frame(data[[i]]))
    
    # load lib
    if (grepl("COLO", toupper(curr_lib))) {
      part1 <- "COREAD"
    }else{
      part1 <- "BRCA"
    }
    part2 <- stringr::str_sub(curr_lib, start = 5, end = 5)
    lib_file_name <- sprintf("%sENCORE_GI_%s_Library_%s.txt", fold_lib, part1, part2)
    
    library[[i]] <- readr::read_tsv(lib_file_name,
                                    col_types = readr::cols(.default = "?", 
                                                            sgRNA1_Chr = "c", 
                                                            sgRNA2_Chr = "c", 
                                                            sgRNA1_WGE_ID = "c", 
                                                            sgRNA2_WGE_ID = "c"), 
                                    show_col_types = FALSE)
    
    # filter library per values in data
    library[[i]] <-  library[[i]][match(data[[i]]$ID, library[[i]]$ID),]
    
    # rename columns
    data[[i]] <- data[[i]] %>%
      dplyr::mutate(lib = curr_lib, .after = ID) %>%
      dplyr::mutate(ID_lib = paste(ID, curr_lib, sep = "_"), .after = lib) %>%
      dplyr::select(-dplyr::ends_with("_logFC"))
    
    tmp <-  library[[i]] %>%
      dplyr::select(ID, sgRNA1_WGE_ID, sgRNA1_WGE_Sequence, sgRNA2_WGE_ID, sgRNA2_WGE_Sequence) %>%
      dplyr::mutate(SEQ_pair = paste0(sgRNA1_WGE_Sequence, "~", sgRNA2_WGE_Sequence))
    
    data[[i]] <- left_join(tmp, data[[i]], by = "ID") %>%
      dplyr::filter(!duplicated(SEQ_pair)) %>%
      dplyr::rename(Note1 = Note, Note2 = MyNote)
    
    # get same order as in data (possible removal of the same SEQ_pair)
    library[[i]] <-  library[[i]][match(data[[i]]$ID, library[[i]]$ID),] %>%
      dplyr::mutate(lib = curr_lib, .after = ID) %>%
      dplyr::mutate(ID_lib = paste(ID, curr_lib, sep = "_"), .after = lib) %>%
      dplyr::mutate(SEQ_pair = paste0(sgRNA1_WGE_Sequence, "~", sgRNA2_WGE_Sequence))
    
    if ("MyNotes" %in% colnames(library[[i]])) {
      library[[i]] <- library[[i]] %>%
        dplyr::select(-MyNotes)
    }
  }
  names(data) <- names(library) <- lib_name
  
  return(list(data = data, library = library))
  
}

# input is a list of output from load_data_rds
get_combined_cneg <- function(res_cneg){
  
  data_cneg <- lapply(res_cneg, function(x) x$data)
  n_cneg <- length(data_cneg)
  data_cneg_logFC <- list()
  for (i in 1:n_cneg) {
    s_tmp <- lapply(data_cneg[[i]], function(x) 
      colnames(x)[grepl("SIDM", colnames(x))])
    data_cneg_logFC[[i]] <- mapply(function(x, y) x[, y], x = data_cneg[[i]], y = s_tmp)
  }
  
  # lib names are the same
  lib_name <- names(data_cneg_logFC[[1]])
  data_cavg <- list()
  
  for (i in 1:length(lib_name)) {
    
    lib_idx <- lib_name[i]
    v_tmp <- lapply(data_cneg_logFC, function(x) x[[lib_idx]])
    
    if (length(unique(sapply(v_tmp, nrow))) > 1) {
      stop("Not the same number of guide pairs!")
    }
    
    list_ID <- lapply(data_cneg, function(x) x[[i]]$ID)
    if (!all(sapply(list_ID, identical, list_ID[[1]]))) {
      stop("Guide pairs do not have the same order!")
    }
    list_colnames <- lapply(v_tmp, colnames)
    if (!all(sapply(list_colnames, identical, list_colnames[[1]]))) {
      stop("Not same samples in cas9negatives")
    }
    data_cavg_logFC <- Reduce('+', v_tmp)/length(v_tmp)
    # use first in the list to annotate
    data_cavg[[i]] <- cbind(data_cneg[[1]][[i]][, !colnames(data_cneg[[1]][[i]]) %in% s_tmp[[i]]], data_cavg_logFC)
  }
  
  names(data_cavg) <- lib_name
  res_cavg <- list(data = data_cavg, library = res_cneg[[1]]$library)
  return(res_cavg)
}

# get sample annotation based on CMP
get_sample_ann <- function(data, file_cmp){
  
  # same for all datasets, use data_cavg
  sample_names <- lapply(data, function(x) 
    colnames(x)[grepl("SID", colnames(x))])
  
  df <- data.frame(model_id_CMP = unlist(sample_names), 
                   lib = unname(unlist(mapply(function(x,y) rep(x,length(y)), 
                                              x = lib_name, 
                                              y = sample_names, SIMPLIFY = T))))
  
  CMP_table <- read_csv(file_cmp, show_col_types = FALSE) %>% 
    dplyr::select(model_id, sample_id, model_name, synonyms, tissue, cancer_type, 
                  tissue_status, COSMIC_ID, BROAD_ID, CCLE_ID) %>%
    dplyr::rename(model_id_CMP = model_id, 
                  sample_id_CMP = sample_id, 
                  model_name_CMP = model_name)
  
  model_encore_table <- dplyr::left_join(df, CMP_table, by = "model_id_CMP") %>%
    mutate(model_name_uppercase = str_replace_all(model_name_CMP, "[-]", "")) %>%
    mutate(model_name_uppercase = toupper(model_name_uppercase))
  model_encore_table <- model_encore_table %>% 
    arrange(lib)
  
  return(model_encore_table)
}

# get common CLs across a list of screens
harmonize_per_CL <- function(list_df, CL_ann){
  
  # rename columns
  mat <- list()
  for (i in 1:length(list_df)) {
    tmp <- as.data.frame(list_df[[i]])
    rownames(tmp) <- tmp$SEQ_pair
    tmp <- tmp[, grepl("SIDM",colnames(tmp))]
    # colnames(tmp) <- str_split_fixed(string = colnames(tmp), pattern = "[_]", n = 4)[,1]
    colnames(tmp) <- CL_ann$model_name_CMP[match(colnames(tmp), CL_ann$model_id_CMP)]
    mat[[i]] <- t(tmp) 
  }
  
  common_CL <- base::Reduce(intersect, lapply(mat, rownames))
  data_common <- lapply(mat, function(x) x[common_CL, ])
  
  names(data_common) <- names(list_df)
  return(data_common)
  
}

# create final tables and save
get_complete_table <- function(list_df, list_matrix){
  
  complete_table <- list()
  df_annot <- mapply(function(x, y) 
    x[match(colnames(y), x$SEQ_pair), !grepl("SIDM", colnames(x))], 
    x = list_df[names(list_matrix)], y = list_matrix, SIMPLIFY = FALSE)
  
  for (i in 1:length(list_matrix)) {
    complete_table[[i]] <- cbind(df_annot[[i]], as.data.frame(t(list_matrix[[i]]))) 
  }
  complete_table <- do.call(rbind, complete_table) %>%
    dplyr::select(-sgRNA1_ID, -sgRNA2_ID)
  
  rownames(complete_table) <- NULL
  return(complete_table)
  
}

ComBatCP <- function(dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                     mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam"),
                     empBayes=TRUE) {
  ## make batch a factor and make a set of indicators for batch
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }
  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  }  ## to be updated soon!
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)  
  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    cat("Using batch =",ref.batch, "as a reference batch (this batch won't change)\n")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")
  
  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  
  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')
  
  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  
  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])
  
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else { 
    B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
  }
  
  ## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }
  
  ## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
    } else {
      var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    }
  } else {
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
    } else {
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
    }
  }
  
  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
  }  
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs){
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  } else{
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE) {
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
    }
  }
  
  if(empBayes){
    ##Find Priors
    gamma.bar <- rowMeans(gamma.hat)
    t2 <- rowVars(gamma.hat)
    a.prior <- apply(delta.hat, 1, sva:::aprior) # FIXME 
    b.prior <- apply(delta.hat, 1, sva:::bprior) # FIXME
    
    ## Plot empirical and parametric priors
    
    if (prior.plots && par.prior) {
      par(mfrow=c(2,2))
      
      ## Top left
      tmp <- density(gamma.hat[1,])
      plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
      xx <- seq(min(tmp$x), max(tmp$x), length=100)
      lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
      
      ## Top Right
      qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
      qqline(gamma.hat[1,], col=2)
      
      ## Bottom Left
      tmp <- density(delta.hat[1,])
      xx <- seq(min(tmp$x), max(tmp$x), length=100)
      tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
      # plot(tmp, typ="l", ylim = c(0, max(tmp$y, tmp1$y)),
      plot(tmp, typ="l",
           main=expression(paste("Density Plot of First Batch ", hat(delta))))
      lines(tmp1, col=2)
      
      ## Bottom Right
      invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
      qqplot(invgam, delta.hat[1,],
             main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
             ylab="Sample Quantiles", xlab="Theoretical Quantiles")
      lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
    }
    
    ## Find EB batch adjustments
    
    gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
    if (par.prior) {
      message("Finding parametric adjustments")
      results <- bplapply(1:n.batch, function(i) {
        if (mean.only) {
          gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
          delta.star <- rep(1, nrow(s.data))
        }
        else {
          temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                               delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                               b.prior[i])
          gamma.star <- temp[1, ]
          delta.star <- temp[2, ]
        }
        list(gamma.star=gamma.star, delta.star=delta.star)
      }, BPPARAM = BPPARAM)
      for (i in 1:n.batch) {
        gamma.star[i,] <- results[[i]]$gamma.star
        delta.star[i,] <- results[[i]]$delta.star
      }
    }
    else {
      message("Finding nonparametric adjustments")
      results <- bplapply(1:n.batch, function(i) {
        if (mean.only) {
          delta.hat[i, ] = 1
        }
        temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                           gamma.hat[i, ], delta.hat[i, ])
        list(gamma.star=temp[1,], delta.star=temp[2,])
      }, BPPARAM = BPPARAM)
      for (i in 1:n.batch) {
        gamma.star[i,] <- results[[i]]$gamma.star
        delta.star[i,] <- results[[i]]$delta.star
      }
    }
  }else{
    #no empirical bayes adjustment:
    gamma.star<-gamma.hat
    delta.star<-delta.hat
  }
  if(!is.null(ref.batch)){
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
  }
  
  ## Normalize the Data ###
  message("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME
  
  ## tiny change still exist when tested on bladder data
  ## total sum of change within each batch around 1e-15 
  ## (could be computational system error).  
  ## Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  
  return(list(correctedData = bayesdata,
              batchDesign = batch.design,
              gamma.star = gamma.star,
              delta.star = delta.star,
              varpool = var.pooled,
              stdmean = stand.mean))
}

# perform combat correction (per guide pairs)
combat_correction <- function(list_df, 
                              CL_ann,
                              save_plot = FALSE,
                              show_plot = TRUE, 
                              save_ext = "png",
                              outfold = NULL){
  
  list_df_h <- harmonize_per_CL(list_df, CL_ann)
  common_pairs <- base::Reduce(intersect, 
                               lapply(list_df, function(x) unique(x$SEQ_pair)))
  
  data_common <- lapply(list_df_h, function(x) t(x[, common_pairs]))
  for (i in 1:length(data_common)) {
    colnames(data_common[[i]]) <- paste(colnames(data_common[[i]]), names(list_df)[i], sep = "_")
  }
  
  tot_mat <- do.call(cbind, data_common)
  ComBat_res <- ComBatCP(
    dat = as.matrix(tot_mat),
    batch = str_split_fixed(colnames(tot_mat), pattern = "_", n = 2)[,2])
  corrected <- as.data.frame(ComBat_res$correctedData)
  
  # annotate combat res and plot
  df_annot <- mapply(function(x, y) 
    x[match(rownames(y), x$SEQ_pair), !grepl("SIDM", colnames(x))], 
    x = list_df, y = data_common, SIMPLIFY = FALSE)
  
  df_annot_unique <- data.frame(SEQ_pair = df_annot[[1]]$SEQ_pair, 
                                Gene_Pair = df_annot[[1]]$Gene_Pair, 
                                Note1 = apply(sapply(df_annot, function(x) x$Note1), 
                                              1, function(y) paste0(sort(unique(y)), collapse = ",")), 
                                Note2 = apply(sapply(df_annot, function(x) x$Note2), 
                                              1, function(y) paste0(sort(unique(y)), collapse = ",")))
  
  colnames(ComBat_res$gamma.star) <- colnames(ComBat_res$delta.star) <- df_annot_unique$SEQ_pair
  rownames(ComBat_res$gamma.star) <- rownames(ComBat_res$delta.star)  <- names(df_annot)
  
  df_ComBat_param_gamma <- melt(ComBat_res$gamma.star, 
                                varnames = c("lib", "SEQ_pair")) %>%
    dplyr::mutate(param = "gamma")
  
  df_ComBat_param_delta <- melt(ComBat_res$delta.star, 
                                varnames = c("lib", "SEQ_pair")) %>%
    dplyr::mutate(param = "delta") 
  df_ComBat_param <- rbind(df_ComBat_param_gamma, df_ComBat_param_delta) %>%
    left_join(df_annot_unique, by = "SEQ_pair")
  
  # plot dist of parameters
  pl_lib <- ggplot(df_ComBat_param, 
                   aes(x = lib, y = value, fill = lib)) + 
    geom_violin() + 
    geom_boxplot(fill = "white", outlier.size = 1, width = 0.2) + 
    theme_bw() + 
    facet_wrap(.~param, scales = "free_y") +
    theme(legend.position = "bottom") +
    xlab("") +
    ylab("ComBat param estimates")
  
  pl_class <- ggplot(df_ComBat_param, 
                     aes(x = Note1, y = value, fill = lib)) + 
    geom_boxplot(outlier.size = 0.5) + 
    theme_bw() + 
    facet_wrap(.~param, scales = "free_x") +
    xlab("") + 
    ylab("ComBat param estimates") +
    theme(legend.position = "bottom") +
    coord_flip()
  
  if (show_plot) {
    print(pl_lib)
    print(pl_class)
  }
  
  if (save_plot) {
    ggsave(filename = sprintf("%sComBat_param_dist_libraries.%s", outfold, save_ext), 
           units = "in", 
           plot = pl_lib, 
           width = 4, 
           height = 4)
    
    ggsave(filename = sprintf("%sComBat_param_dist_GPclass.%s", outfold, save_ext), 
           units = "in", 
           plot = pl_class, 
           width = 6, 
           height = 5)
  }
  
  return(list(raw = tot_mat,
              corrected = corrected,
              ComBat_res = ComBat_res))
}


# PC plot for common pairs
# plot before and after combat correction
pca_commonpairs_function <- function(list_df, 
                                     CL_ann,
                                     save_plot = FALSE,
                                     save_ext = "png",
                                     show_plot = TRUE, 
                                     outfold = NULL){
  
  res <- combat_correction(list_df, 
                           CL_ann,
                           show_plot = FALSE)
  
  corrected_common <- res$corrected
  raw_common <- res$raw
  common_pairs <- rownames(raw_common)
  
  # pc from corrected data
  pca_common <- prcomp(t(corrected_common), scale. = TRUE)
  pc_corr <- data.frame(pca_common$x,
                        CL = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,1],
                        lib = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,2])
  pc_corr$lib <- factor(pc_corr$lib)
  pc_corr$CL_label <- ""
  pc_corr$CL_label[grepl("1", pc_corr$lib)] <- pc_corr$CL[grepl("1", pc_corr$lib)] 
  pc_corr$CL <- factor(pc_corr$CL)
  
  pl1 <- ggplot(pc_corr, aes(x = PC1,
                             y = PC2,
                             col = CL,
                             label = CL_label,
                             group = CL)) +
    geom_point(size = 2, aes(shape = lib)) +
    geom_polygon(aes(fill = CL), alpha = 0.2, size = 0.1) +
    # geom_line(alpha = 0.5) + 
    geom_text_repel(size = 3, max.overlaps = Inf) +
    guides(color = "none", fill = "none") +
    theme_bw() +
    ggtitle("ComBat corrected",
            subtitle = paste0("N. common guide pairs: ", length(common_pairs)))
  
  # pc from raw data
  pca_common <- prcomp(t(raw_common), scale. = TRUE)
  pc_raw <- data.frame(pca_common$x,
                       CL = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,1],
                       lib = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,2])
  pc_raw$lib <- factor(pc_raw$lib)
  pc_raw$CL_label <- ""
  pc_raw$CL_label[grepl("1", pc_raw$lib)] <- pc_raw$CL[grepl("1", pc_raw$lib)] 
  pc_raw$CL <- factor(pc_raw$CL)
  
  pl2 <- ggplot(pc_raw, aes(x = PC1,
                            y = PC2,
                            col = CL,
                            label = CL_label,
                            group = CL)) +
    geom_point(size = 2, aes(shape = lib)) +
    geom_polygon(aes(fill = CL), alpha = 0.2, size = 0.1) +
    # geom_line(alpha = 0.5) + 
    geom_text_repel(size = 3, max.overlaps = Inf) +
    guides(color = "none", fill = "none") +
    theme_bw() +
    ggtitle("Raw logFC",
            subtitle = paste0("N. common guide pairs: ", length(common_pairs)))
  
  pl <- ggpubr::ggarrange(plotlist = list(pl2, pl1), ncol = 2, common.legend = TRUE)
  
  if (show_plot) {
    print(pl)
  }
  
  if (save_plot) {
    ggsave(filename = sprintf("%sPC1_2_raw_and_corrected.%s", outfold, save_ext), 
           units = "in", 
           plot = pl, 
           width = 8, 
           height = 5.5)
  }
  
  return(list(common_pairs = common_pairs,
              pc_corr = pc_corr,
              pc_raw = pc_raw))
  
}

# get median distance
dist_commonpairs <- function(mat_common){
  
  dist_mat <- as.matrix(dist(t(mat_common),method = 'euclidean'))
  # cor_res <- cor(mat_common, method = "pearson")
  libs <- str_split_fixed(colnames(dist_mat), pattern = "_", n = 2)[,2]
  libs_name <- unique(libs)
  CLs <- str_split_fixed(colnames(dist_mat), pattern = "_", n = 2)[,1]
  CLs_names <- unique(CLs)
  
  inner_libs <- sapply(libs_name, function(x) median(dist_mat[libs == x, libs == x][upper.tri(dist_mat[libs == x, libs == x])]))
  outer_libs <- sapply(libs_name, function(x) median(dist_mat[libs == x, libs != x]))
  
  inner_CLs <- sapply(CLs_names, function(x) median(dist_mat[CLs == x, CLs == x][upper.tri(dist_mat[CLs == x, CLs == x])]))
  outer_CLs <- sapply(CLs_names, function(x) median(dist_mat[CLs == x, CLs != x]))
  
  return(list(libs = data.frame(name = names(inner_libs), inner = inner_libs, outer = outer_libs), 
              CLs = data.frame(name = names(inner_CLs), inner = inner_CLs, outer = outer_CLs)))
}

# get ordered list
sorted_dist <- function(mat_common){
  
  dist_matrix <- as.matrix(dist(t(mat_common),
                                method = 'euclidean'))
  
  upper_tri_indices <- which(upper.tri(dist_matrix, diag = FALSE), arr.ind = TRUE)
  
  df_dist <- data.frame(
    sample1 = rownames(dist_matrix)[upper_tri_indices[, 1]],
    sample2 = colnames(dist_matrix)[upper_tri_indices[, 2]],
    dist = dist_matrix[upper_tri_indices]
  ) %>% 
    dplyr::arrange(dist) %>%
    dplyr::mutate(sample1_CL = str_split_fixed(sample1, pattern = "_", n = 2)[,1], 
                  sample2_CL = str_split_fixed(sample2, pattern = "_", n = 2)[,1],
                  sample1_lib = str_split_fixed(sample1, pattern = "_", n = 2)[,2], 
                  sample2_lib = str_split_fixed(sample2, pattern = "_", n = 2)[,2]) %>%
    dplyr::mutate(same_lib = sample1_lib == sample2_lib, 
                  same_CL = sample1_CL == sample2_CL)
  
  return(df_dist)  
}


plot_dist_PPV <- function(list_df,  
                          CL_ann,
                          save_plot = FALSE,
                          save_ext = "png",
                          outfold = NULL){
  
  res_combat <- combat_correction(list_df, 
                                  CL_ann,
                                  save_plot = FALSE,
                                  show_plot = FALSE, 
                                  outfold = NULL)
  
  dist_sorted_adj <- sorted_dist(mat_common = res_combat$corrected)
  df_PPV_adj <- data.frame(
    CL = cumsum(dist_sorted_adj$same_CL)/1:nrow(dist_sorted_adj), 
    lib = cumsum(dist_sorted_adj$same_lib)/1:nrow(dist_sorted_adj),
    id = 1:nrow(dist_sorted_adj)
  )
  df_PPV_adj <- melt(data = df_PPV_adj, 
                     id.vars = "id",
                     variable.name = "type", 
                     value.name = "PPV")
  
  dist_sorted_raw <- sorted_dist(mat_common = res_combat$raw)
  df_PPV_raw <- data.frame(
    CL = cumsum(dist_sorted_raw$same_CL)/1:nrow(dist_sorted_raw), 
    lib = cumsum(dist_sorted_raw$same_lib)/1:nrow(dist_sorted_raw),
    id = 1:nrow(dist_sorted_adj)
  )
  df_PPV_raw <- melt(data = df_PPV_raw, 
                     id.vars = "id",
                     variable.name = "type", 
                     value.name = "PPV")
  
  pl1 <- ggplot(df_PPV_raw, aes(x = id,
                                y = PPV,
                                color = type)) +
    geom_line() +
    xlab("K closest sample pair") + 
    ylab("PPV = n. same CL or lib / K") + 
    theme_bw() + 
    theme(legend.title = element_blank()) +
    ggtitle("Raw")
  
  pl2 <- ggplot(df_PPV_adj, aes(x = id,
                                y = PPV,
                                color = type)) +
    geom_line() +
    xlab("K closest sample pair") + 
    ylab("PPV = n. same CL or lib / K") + 
    theme_bw() + 
    theme(legend.title = element_blank()) +
    ggtitle("ComBat Corrected")
  
  pl <- ggpubr::ggarrange(plotlist = list(pl1, pl2), ncol = 2, 
                          common.legend = TRUE)
  print(pl)
  if (save_plot) {
    ggsave(filename = sprintf("%sPPV_CL-lib.%s", outfold, save_ext), 
           units = "in", 
           plot = pl, 
           width = 7, 
           height = 4)
  }
}


plot_dist_commonpairs <- function(list_df,
                                  CL_ann,
                                  save_plot = FALSE, 
                                  save_ext = "png",
                                  outfold = NULL){
  
  res_combat <- combat_correction(list_df, 
                                  CL_ann,
                                  save_plot = FALSE,
                                  show_plot = FALSE, 
                                  outfold = NULL)
  
  dist_raw <- dist_commonpairs(res_combat$raw)
  dist_combat <- dist_commonpairs(res_combat$corrected)
  
  CL_dist <- data.frame(name = rep(dist_raw$CLs$name, 2), dist = c(dist_raw$CLs$inner, dist_raw$CLs$outer), 
                        type_dist = c(rep("same", nrow(dist_raw$CLs)), rep("different", nrow(dist_raw$CLs))), 
                        type_combat = "Raw")
  CL_dist <- rbind(CL_dist, data.frame(name = rep(dist_combat$CLs$name, 2), dist = c(dist_combat$CLs$inner, dist_combat$CLs$outer), 
                                       type_dist = c(rep("same", nrow(dist_combat$CLs)), rep("different", nrow(dist_combat$CLs))), 
                                       type_combat = "ComBat Corrected"))
  CL_dist$type_combat <- factor(CL_dist$type_combat, levels = c("Raw", "ComBat Corrected"))
  
  pl1 <- ggplot(CL_dist , aes(x = type_combat,
                              y = dist,
                              fill = type_dist)) +
    geom_boxplot(alpha = 0.4) +
    geom_jitter(position = position_dodge(width = 0.7)) +
    xlab("") + 
    ylab("Median Eucl. distance") + 
    theme_bw() + 
    ggtitle("CLs")
  
  libs_dist <- data.frame(name = rep(dist_raw$libs$name, 2), dist = c(dist_raw$libs$inner, dist_raw$libs$outer), 
                          type_dist = c(rep("same", nrow(dist_raw$libs)), rep("different", nrow(dist_raw$libs))), 
                          type_combat = "Raw")
  libs_dist <- rbind(libs_dist, data.frame(name = rep(dist_combat$libs$name, 2), dist = c(dist_combat$libs$inner, dist_combat$libs$outer), 
                                           type_dist = c(rep("same", nrow(dist_combat$libs)), rep("different", nrow(dist_combat$libs))), 
                                           type_combat = "ComBat Corrected"))
  libs_dist$type_combat <- factor(libs_dist$type_combat, levels = c("Raw", "ComBat Corrected"))
  pl2 <- ggplot(libs_dist , aes(x = type_combat,
                                y = dist,
                                fill = type_dist)) +
    geom_boxplot(alpha = 0.4) +
    geom_jitter(position = position_dodge(width = 0.7)) +
    xlab("") + 
    ylab("Median Eucl. distance") + 
    theme_bw() + 
    ggtitle("Libraries")
  pl <- ggpubr::ggarrange(plotlist = list(pl2, pl1), ncol = 2, common.legend = TRUE)
  print(pl)
  
  if (save_plot) {
    ggsave(filename = sprintf("%sinside_and_outside_distance_same_CL-lib.%s", outfold, save_ext), 
           units = "in", 
           plot = pl, 
           width = 7, 
           height = 4)
  }
  
  return(list(CL = CL_dist, 
              lib = libs_dist))
  
}

# can we estimate gamma and delta from the closest elements?
# validation based on cohen's d inside and outside kNN
get_stat_closest <- function(id_lib, 
                             matrix_data, 
                             combat_param, 
                             kNN){
  
  # explain why euclidean distance makes sense
  dist_guides <- as.matrix(dist(t(matrix_data[[id_lib]])))
  # for each guide, find closest guides
  closest_guides <- apply(dist_guides, 1, function(x) order(x)[1:(kNN + 1)], simplify = FALSE)
  closest_guides_dist <- apply(dist_guides, 1, function(x) sort(x)[1:(kNN + 1)], simplify = FALSE)
  
  gamma_closest <- lapply(closest_guides, 
                          function(x) combat_param$gamma.star[id_lib,x])
  gamma_compl <- lapply(closest_guides, 
                        function(x) combat_param$gamma.star[id_lib,-x])
  delta_closest <- lapply(closest_guides, 
                          function(x) combat_param$delta.star[id_lib,x])
  delta_compl <- lapply(closest_guides, 
                        function(x) combat_param$delta.star[id_lib,-x])
  
  ###
  n_pairs <- ncol(matrix_data[[id_lib]])
  cohensd_gamma <- sapply(1:n_pairs, function(x) 
    tryCatch(cohens_d(
      x = abs(combat_param$gamma.star[id_lib,x] - gamma_closest[[x]]), 
      y = abs(combat_param$gamma.star[id_lib,x] - gamma_compl[[x]]), 
      pooled_sd = FALSE)$Cohens_d, 
      error = function(e){NA}))
  
  ttest_gamma <- sapply(1:n_pairs, function(x) 
    tryCatch(t.test(
      x = abs(combat_param$gamma.star[id_lib,x] - gamma_closest[[x]]),
      y = abs(combat_param$gamma.star[id_lib,x] - gamma_compl[[x]]), 
      alternative = "less")$p.value, 
      error = function(e){NA}))
  
  cohensd_delta <- sapply(1:n_pairs, function(x) 
    tryCatch(cohens_d(
      x = abs(combat_param$delta.star[id_lib,x] - delta_closest[[x]]), 
      y = abs(combat_param$delta.star[id_lib,x] - delta_compl[[x]]), 
      pooled_sd = FALSE)$Cohens_d, 
      error = function(e){NA}))
  
  ttest_delta <- sapply(1:n_pairs, function(x) 
    tryCatch(t.test(
      x = abs(combat_param$delta.star[id_lib,x] - delta_closest[[x]]), 
      y = abs(combat_param$delta.star[id_lib,x] - delta_compl[[x]]), 
      alternative = "less")$p.value,
      error = function(e){NA}))
  
  df_out <- data.frame(SEQ_pair = rep(colnames(matrix_data[[id_lib]]), 2), 
                       cohens_d = c(cohensd_gamma, cohensd_delta), 
                       ttest_pval = c(ttest_gamma, ttest_delta), 
                       type = c(rep("gamma", n_pairs), rep("delta", n_pairs)), 
                       id_lib = id_lib)
  
  df_out$kNN <- kNN
  
  return(df_out)
  
}

# can we estimate gamma and delta from the closest elements?
# validation based on auc using ED distance between param as predictor
compute_auc <- function(id_lib, 
                        matrix_data, 
                        combat_param, 
                        kNN){
  
  # distance of logFCs across CLs
  dist_guides <- as.matrix(dist(t(matrix_data[[id_lib]])))
  closest_guides <- apply(dist_guides, 1, function(x) order(x)[1:(kNN + 1)], simplify = FALSE)
  n_pairs <- ncol(matrix_data[[id_lib]])
  
  # distance of the gamma,delta parameters across guide pairs 
  param_vector <- data.frame(gamma = combat_param$gamma.star[id_lib,], 
                             delta = combat_param$delta.star[id_lib,])
  dist_param <- as.matrix(dist(param_vector, method = "euclidean"))
  
  # response = kNN of the guide pair
  # predictor = distance between the guide pair and all the other guide pairs
  auc_res <- sapply(1:n_pairs, function(x)
    as.numeric(auc(roc(
      response = as.numeric(1:n_pairs %in% closest_guides[[x]]), 
      predictor = unname(dist_param[,x]), 
      quiet = TRUE, 
      direction = ">"))))
  
  df_out <- data.frame(SEQ_pair = colnames(matrix_data[[id_lib]]),
                       auc = auc_res, 
                       id_lib = id_lib, 
                       kNN = kNN)
  return(df_out)
}

# combine all validation strategies:
validate_NN_approximation <- function(list_df, 
                                      CL_ann,
                                      save_plot = FALSE,
                                      save_ext = "png",
                                      outfold = NULL){
  
  res <- combat_correction(list_df, 
                           CL_ann = CL_ann,
                           save_plot = FALSE, 
                           show_plot = TRUE)
  
  corrected_common <- res$corrected
  raw_common <- res$raw
  common_pairs <- rownames(raw_common)
  matrix_data <- lapply(harmonize_per_CL(list_df, CL_ann), function(x) 
    x[,common_pairs])
  
  # varying kNN
  kNN_list <- c(5, seq(10, 50, by = 10))
  df_kNN <- data.frame()
  # Validation 1: cohens'D  between neightbours of gp i vs not neighb 
  # (based on euclidean distance)
  for (kNN_p in kNN_list) {
    print(kNN_p)
    tmp <- lapply(1:length(list_df), function(x)
      get_stat_closest(id_lib = x, 
                       matrix_data = matrix_data, 
                       combat_param = res$ComBat_res, 
                       kNN = kNN_p))  
    df_kNN <- rbind(do.call(rbind, tmp), df_kNN)
  }
  
  pl <- ggplot(df_kNN, aes(x = as.factor(kNN),
                           y = cohens_d)) +
    geom_boxplot() +
    facet_wrap(.~type) +
    xlab("k Nearest Neigh.") +
    ylab("Cohen's d:\n |param_i - neigh.| VS |param_i - not neigh.|") +
    theme_bw()
  print(pl)
  if (save_plot) {
    ggsave(plot = pl, 
           filename = sprintf("%svalidation_kNN_cohensd.%s", outfold, save_ext), 
           width = 6, 
           height = 4)    
  }
  
  ### Validation 2:
  # create TPR, FPR. Compare the selection of closest kNN with the "closest" param
  auc_kNN <- data.frame()
  for (kNN_p in kNN_list) {
    print(kNN_p)
    tmp <- lapply(1:length(list_df), function(x)
      compute_auc(id_lib = x, 
                  matrix_data = matrix_data, 
                  combat_param = res$ComBat_res, 
                  kNN = kNN_p))  
    auc_kNN <- rbind(do.call(rbind, tmp), auc_kNN)
  }
  
  pl <- ggplot(auc_kNN, aes(x = as.factor(kNN),
                            y = auc)) +
    geom_boxplot() +
    xlab("k Nearest Neigh.") +
    ylab("AUROC") +
    theme_bw()
  print(pl)
  if (save_plot) {
    ggsave(plot = pl, 
           filename = sprintf("%svalidation_kNN_AUC.%s", outfold, save_ext), 
           width = 4, 
           height = 4)
  }
  
  
  return(list(cohensd_kNN = df_kNN, 
              auc_kNN = auc_kNN))
  
}

# approximate combat param via kNN
param_approx_kNN <- function(
    matrix_data_list, 
    gamma_matrix, 
    delta_matrix, 
    kNN, 
    outfold = NULL, 
    save_plot = FALSE, 
    save_ext = "png",
    show_plot = TRUE, 
    list_df = NULL){
  
  res <- list()
  for (id in 1:length(matrix_data_list)) {
    
    matrix_data <- matrix_data_list[[id]]
    gamma_param <- gamma_matrix[id,]
    delta_param <- delta_matrix[id,]
    
    common_pairs <- names(gamma_param)
    unique_pairs <- setdiff(colnames(matrix_data), common_pairs)
    
    # distance between common pairs and unique pairs only
    dist_logFCs <- sqrt(dist2(t(matrix_data[,common_pairs]),
                              t(matrix_data[,unique_pairs])))
    
    closest_guides <- apply(dist_logFCs, 2, function(x) 
      order(x)[1:kNN], 
      simplify = TRUE)
    
    gamma_approx_unique <- apply(closest_guides, 2, function(x)
      mean(gamma_param[common_pairs[x]])
    )
    
    delta_approx_unique <- apply(closest_guides, 2, function(x)
      mean(delta_param[common_pairs[x]])
    )
    res[[id]] <- list(gamma = gamma_approx_unique, 
                      delta = delta_approx_unique)
  }
  
  
  if (show_plot) {
    
    param_approx <- res
    libs <- names(matrix_data_list)
    
    df_annot <- mapply(function(x, y) 
      y[match(names(x$delta), y$SEQ_pair), !grepl("SIDM", colnames(y))], 
      x = param_approx, y = list_df, SIMPLIFY = FALSE)
    df_annot <- do.call(rbind, df_annot)
    
    gamma_tmp <- lapply(param_approx, function(x) x$gamma)
    df_pl_gamma <- data.frame(value = unlist(gamma_tmp), 
                              param = "gamma", 
                              df_annot)
    
    delta_tmp <- lapply(param_approx, function(x) x$delta)
    df_pl_delta <- data.frame(value = unlist(delta_tmp), 
                              param = "delta", 
                              df_annot)
    df_pl <- rbind(df_pl_delta, df_pl_gamma)
    
    # plot dist of parameters
    pl_lib <- ggplot(df_pl, 
                     aes(x = lib, y = value, fill = lib)) + 
      geom_violin() + 
      geom_boxplot(fill = "white", outlier.size = 1, width = 0.2) + 
      theme_bw() + 
      facet_wrap(.~param, scales = "free_y") +
      theme(legend.position = "bottom") +
      xlab("") +
      ylab("Param KNN approximation")
    
    pl_class <- ggplot(df_pl, 
                       aes(x = Note1, y = value, fill = lib)) + 
      geom_boxplot(outlier.size = 0.5) + 
      theme_bw() + 
      facet_wrap(.~param, scales = "free_x") +
      xlab("") + 
      ylab("Param KNN approximation") +
      theme(legend.position = "bottom") +
      coord_flip()
    
    print(pl_lib)
    print(pl_class)
    
    if (save_plot) {
      ggsave(filename = sprintf("%skNNApprox_param_dist_libraries.%s", outfold, save_ext), 
             units = "in", 
             plot = pl_lib, 
             width = 4, 
             height = 4)
      
      ggsave(filename = sprintf("%skNNApprox_param_dist_GPclass.%s", outfold, save_ext), 
             units = "in", 
             plot = pl_class, 
             width = 6, 
             height = 5)
    }
  }
  
  return(res)
  
}

# distance between each element in a matrix
dist2 <- function(X,C) {
  # from SNFtools
  
  ndata <- nrow(X)
  ncentres <- nrow(C)
  
  sumsqX <- rowSums(X^2)
  sumsqC <- rowSums(C^2)
  
  XC <- 2 * (X %*% t(C))
  res <- matrix(rep(sumsqX, times = ncentres), ndata, ncentres) + 
    t(matrix(rep(sumsqC, times = ndata), ncentres, ndata)) - XC
  res[res < 0] <- 0
  
  return(res)
}


# adjust data, combat estimates per guide pairs
adjust_alldata_kNN <- function(list_df, 
                               CL_ann = CL_ann, 
                               kNN, 
                               outfold = NULL, 
                               save_plot = FALSE, 
                               save_ext = "png",
                               show_plot = TRUE){
  
  combat_res_all <- combat_correction(list_df = list_df, 
                                      CL_ann = CL_ann,
                                      outfold = outfold, 
                                      save_plot = save_plot, 
                                      save_ext = save_ext,
                                      show_plot = show_plot)
  
  combat_res_all$common_pairs <- rownames(combat_res_all$ComBat_res$correctedData)
  
  combat_res <- combat_res_all$ComBat_res
  list_logFC <- harmonize_per_CL(list_df, CL_ann)
  n_guide_common <- nrow(combat_res$correctedData)
  
  # get combat param
  batch.design <- combat_res$batchDesign
  gamma.star <- combat_res$gamma.star # additive batch effect
  delta.star <- combat_res$delta.star # multiplicative batch effect
  stand.mean <- combat_res$stdmean
  var.pooled <- combat_res$varpool
  
  # standardize each guide pair, Mean and Sd computed across CLs
  stand.mean_all <- lapply(list_logFC, function(x) rowMeans(t(x)) %*% t(rep(1, ncol(t(x)))))
  var.pooled_all <- lapply(list_logFC, function(x) rowSds(t(x)) %*% t(rep(1, ncol(t(x)))))
  # NOTE:
  # we cannot use stand.mean and var.pooled (stdPrior) because we need to estimate guide pairs not in common
  # but correlation is high for common pairs with this and ComBatCP strategy 
  # >0.95 for stand.mean and >0.84 for var.pooled
  s.data_all <- mapply(function(x,y,z) (t(x) - y)/z, 
                       x = list_logFC, 
                       y = stand.mean_all, 
                       z = var.pooled_all, 
                       SIMPLIFY = FALSE)
  
  # correct based on bayes estimates
  mean_common <- mapply(function(x, y) 
    matrix(gamma.star[x,], 
           nrow = n_guide_common, 
           ncol = nrow(y)), 
    x = 1:nrow(gamma.star), y = list_logFC, SIMPLIFY = FALSE)
  
  var_common <- mapply(function(x, y) 
    sqrt(delta.star[x,]) %*% t(rep(1,nrow(y))), 
    x = 1:nrow(delta.star), y = list_logFC, 
    SIMPLIFY = FALSE)
  
  # for those not in common, use kNN
  param_approx <- param_approx_kNN(matrix_data_list = list_logFC, 
                                   gamma_matrix = gamma.star,
                                   delta_matrix = delta.star, 
                                   kNN = kNN,
                                   outfold = outfold, 
                                   save_ext = save_ext,
                                   save_plot = save_plot, 
                                   show_plot = show_plot, 
                                   list_df = list_df)
  
  mean_unique <- var_unique <- list()
  mean_all <- var_all <- list()
  
  for (id in 1:length(list_logFC)) {
    
    mean_unique[[id]] <- matrix(param_approx[[id]]$gamma, 
                                nrow = length(param_approx[[id]]$gamma), 
                                ncol = nrow(list_logFC[[id]]))
    
    var_unique[[id]] <- matrix(param_approx[[id]]$delta, 
                               nrow = length(param_approx[[id]]$delta), 
                               ncol = nrow(list_logFC[[id]]))
    
    rownames(mean_common[[id]]) <- colnames(gamma.star)
    rownames(var_common[[id]]) <- colnames(delta.star)
    rownames(mean_unique[[id]]) <- names(param_approx[[id]]$gamma)
    rownames(var_unique[[id]]) <- names(param_approx[[id]]$delta)
    
    # put back together
    var_all[[id]] <- rbind(var_common[[id]], var_unique[[id]])
    mean_all[[id]] <- rbind(mean_common[[id]], mean_unique[[id]])
    # order properly
    var_all[[id]] <- var_all[[id]][rownames(s.data_all[[id]]),]
    mean_all[[id]] <- mean_all[[id]][rownames(s.data_all[[id]]),]
    
  }
  
  
  adjusted.data_all <- mapply(function(x,y,z) (x - y)/z, 
                              x = s.data_all, 
                              y =  mean_all, 
                              z = var_all, 
                              SIMPLIFY = FALSE)
  
  # add original mean and variance back
  adjusted.data_all <- mapply(function(x,y,z) t( (x*z) + y ), 
                              x = adjusted.data_all, 
                              y = stand.mean_all, 
                              z = var.pooled_all, 
                              SIMPLIFY = FALSE)
  
  # save param
  gamma_mean <- lapply(mean_all, function(x) x[,1])
  names(gamma_mean) <- names(list_df)
  delta_var <- lapply(var_all, function(x) x[,1])
  names(delta_var) <- names(list_df)
  
  return(list(adj = adjusted.data_all, 
              original = list_logFC, 
              combat = combat_res_all, 
              param_notmatching = param_approx, 
              param_all = list(gamma_mean = gamma_mean, delta_var = delta_var)))
  
}

# plot distributions per CL
plot_CL_distribution <- function(original, 
                                 adjusted, 
                                 list_df,
                                 common_pairs, 
                                 outfold = NULL, 
                                 show_plot = TRUE,
                                 save_plot = FALSE, 
                                 save_ext = "png"
){
  
  CL_names <- rownames(original[[1]])
  # get annotation
  df_annot <- mapply(function(x, y) 
    x[match(colnames(y), x$SEQ_pair), !grepl("SIDM", colnames(x))], 
    x = list_df[names(original)], y = original, SIMPLIFY = FALSE)
  
  df_annot_common <- lapply(df_annot, function(x) x[match(common_pairs, x$SEQ_pair),])
  # get unique
  df_annot_common_unique <- data.frame(SEQ_pair = df_annot_common[[1]]$SEQ_pair, 
                                       Gene_Pair = df_annot_common[[1]]$Gene_Pair, 
                                       Note1 = apply(sapply(df_annot_common, function(x) x$Note1), 
                                                     1, function(y) paste0(sort(unique(y)), collapse = ",")), 
                                       Note2 = apply(sapply(df_annot_common, function(x) x$Note2), 
                                                     1, function(y) paste0(sort(unique(y)), collapse = ",")))
  ## plot everything  
  pl <- list()
  df_CL <- list()
  for (i in 1:length(CL_names)) {
    
    CL_name <- CL_names[i]
    print(sprintf("%s", CL_name))
    
    df_tmp <- mapply(function(x,y,z) 
      cbind(z, data.frame(logFC_or = x[CL_name,, drop = T], logFC_adj = y[CL_name,, drop = T])), 
      x = original, y = adjusted, z = df_annot, SIMPLIFY = F)
    
    df_CL[[i]] <- bind_rows(df_tmp)
    print(dim(df_CL[[i]]))
    
    df_CL[[i]] <- left_join(suffix = c("lib_spec", "lib_common"), 
                            df_CL[[i]], 
                            df_annot_common_unique, 
                            by = "SEQ_pair")
    df_CL[[i]]$Note1 <- df_CL[[i]]$Note1lib_common
    df_CL[[i]]$Note1[is.na(df_CL[[i]]$Note1lib_common)] <- df_CL[[i]]$Note1lib_spec[is.na(df_CL[[i]]$Note1lib_common)]
    
    pl1 <- ggplot(df_CL[[i]], 
                  aes(x = Note1, y = logFC_adj, fill = lib)) + 
      geom_boxplot(outlier.size = 0.5) + 
      theme_bw() + 
      xlab("") +
      coord_flip() + 
      theme(axis.text.y = element_blank()) + 
      ggtitle(sprintf("%s: Corrected", CL_name))
    
    pl2 <- ggplot(df_CL[[i]],
                  aes(x = Note1, y = logFC_or, fill = lib)) + 
      geom_boxplot(outlier.size = 0.5) + 
      theme_bw() + 
      xlab("") +
      coord_flip() + 
      ggtitle(sprintf("%s: Raw", CL_name))
    
    pl3 <- ggplot(subset(df_CL[[i]], SEQ_pair %in% common_pairs),
                  aes(x = Note1, y = logFC_adj, fill = lib)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_bw() +
      xlab("") +
      coord_flip() +
      theme(axis.text.y = element_blank()) + 
      ggtitle(sprintf("%s: Corrected (common pairs)", CL_name))
    
    pl4 <- ggplot(subset(df_CL[[i]], SEQ_pair %in% common_pairs),
                  aes(x = Note1, y = logFC_or, fill = lib)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_bw() +
      xlab("") +
      coord_flip() +
      ggtitle(sprintf("%s: Raw (common pairs)", CL_name))
    
    pl[[i]] <- ggpubr::ggarrange(plotlist = list(pl4, pl3, pl2, pl1), 
                                 heights = c(0.8,1),
                                 widths = c(1, 0.6),
                                 ncol = 2, 
                                 nrow = 2, 
                                 common.legend = TRUE)
    if (show_plot) {
      print(pl[[i]])
    }
    
  }
  
  # plot dist across all CLs
  df_tot <- do.call(rbind, df_CL)  
  pl1 <- ggplot(df_tot, 
                aes(x = Note1, y = logFC_adj, fill = lib)) + 
    geom_boxplot(outlier.size = 0.2) + 
    theme_bw() + 
    xlab("") +
    ylab("logFC") +
    coord_flip() + 
    theme(axis.text.y = element_blank()) + 
    ggtitle(sprintf("%s: Corrected", "All CLs"))
  
  pl2 <- ggplot(df_tot,
                aes(x = Note1, y = logFC_or, fill = lib)) + 
    geom_boxplot(outlier.size = 0.2) + 
    theme_bw() + 
    xlab("") +
    ylab("logFC") +
    coord_flip() + 
    theme(axis.text.y = element_text(size = 10)) + 
    ggtitle(sprintf("%s: Raw", "All CLs"))
  
  pl3 <- ggplot(subset(df_tot, SEQ_pair %in% common_pairs), 
                aes(x = Note1lib_common, y = logFC_adj, fill = lib)) + 
    geom_boxplot(outlier.size = 0.2) + 
    theme_bw() + 
    xlab("") +
    ylab("logFC") +
    coord_flip() + 
    theme(axis.text.y = element_blank()) + 
    ggtitle(sprintf("%s: Corrected (common pairs)", "All CLs"))
  
  pl4 <- ggplot(subset(df_tot, SEQ_pair %in% common_pairs),
                aes(x = Note1lib_common, y = logFC_or, fill = lib)) + 
    geom_boxplot(outlier.size = 0.2) + 
    theme_bw() + 
    xlab("") +
    ylab("logFC") + 
    coord_flip() + 
    theme(axis.text.y = element_text(size = 10)) + 
    ggtitle(sprintf("%s: Raw (common pairs)", "All CLs"))
  
  pl_class <- ggpubr::ggarrange(plotlist = list(pl4, pl3, pl2, pl1), 
                                heights = c(0.8,1),
                                widths = c(1, 0.6),
                                ncol = 2, 
                                nrow = 2, 
                                align = "h", 
                                common.legend = TRUE)
  
  if (save_plot) {
    ggsave(filename = sprintf("%sALLCLs_distr_perClass.%s", outfold, save_ext), 
           units = "in", 
           plot = pl_class, 
           width = 8, 
           height = 8)
  }
  
  pl1 <- ggplot(df_tot, 
                aes(x = lib, y = logFC_adj, fill = lib)) + 
    geom_boxplot(outlier.size = 0.5) + 
    theme_bw() + 
    xlab("") +
    ylab("logFC") +
    ggtitle(sprintf("%s: Corrected", "All CLs"))
  
  pl2 <- ggplot(df_tot,
                aes(x = lib, y = logFC_or, fill = lib)) + 
    geom_boxplot(outlier.size = 0.5) + 
    theme_bw() + 
    xlab("") +
    ylab("logFC") +
    ggtitle(sprintf("%s: Raw", "All CLs"))
  
  pl3 <- ggplot(subset(df_tot, SEQ_pair %in% common_pairs), 
                aes(x = lib, y = logFC_adj, fill = lib)) + 
    geom_boxplot(outlier.size = 0.5) + 
    theme_bw() + 
    xlab("") +
    ylab("logFC") +
    ggtitle(sprintf("%s: Corrected (common pairs)", "All CLs"))
  
  pl4 <- ggplot(subset(df_tot, SEQ_pair %in% common_pairs),
                aes(x = lib, y = logFC_or, fill = lib)) + 
    geom_boxplot(outlier.size = 0.5) + 
    theme_bw() + 
    xlab("") +
    ylab("logFC") +
    ggtitle(sprintf("%s: Raw (common pairs)", "All CLs"))
  
  pl_lib <- ggpubr::ggarrange(plotlist = list(pl4, pl3, pl2, pl1),
                              ncol = 2, 
                              nrow = 2, 
                              align = "h", 
                              common.legend = TRUE)
  
  
  if (show_plot) {
    print(pl_lib)
  }
  
  if (save_plot) {
    ggsave(filename = sprintf("%sALLCLs_distr_perLib.%s", outfold, save_ext), 
           units = "in", 
           plot = pl_lib, 
           width = 6, 
           height = 7)
  }
  
  return(df_tot)
}

test_distributions_per_class_allCLs <- function(df, 
                                                outfold = NULL, 
                                                save_plot = FALSE, 
                                                save_ext = "png"){
  
  df <- df %>%
    dplyr::filter(!Note1 %in% "AnchorCombinations")
  class_names <- unique(df$Note1)
  
  ktest_adj <- sapply(class_names, function(x) kruskal.test(x = df$logFC_adj[df$Note1 %in% x], 
                                                            g = df$lib[df$Note1 %in% x])$p.value)
  ktest_adj <- data.frame(pvalue = unname(ktest_adj), class = names(ktest_adj), type = "ComBat Adjusted")
  
  ktest_or <- sapply(class_names, function(x) kruskal.test(x = df$logFC_or[df$Note1 %in% x], 
                                                           g = df$lib[df$Note1 %in% x])$p.value)
  ktest_or <- data.frame(pvalue = unname(ktest_or), class = names(ktest_or), type = "Raw")
  
  ktest_all <- rbind(ktest_adj, ktest_or) %>% 
    mutate(log10p = -log10(pvalue + .Machine$double.xmin))
  
  pl <- ggplot(ktest_all, aes(x = class,
                              y = log10p,
                              col = type)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "grey50") +
    # geom_bar(stat = "identity", position = position_dodge()) + 
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    # geom_text_repel(size = 3, max.overlaps = Inf) +
    theme_bw() + 
    theme(legend.position = "bottom", legend.title = element_blank()) + 
    xlab("") + 
    ylab("-log10 P-value Kruskal-Wallis test") + 
    ggtitle("All CLs") + 
    coord_flip()
  print(pl)
  
  
  if (save_plot) {
    ggsave(filename = sprintf("%skruskaltest_perclass_ALLCLs.%s", outfold, save_ext), 
           units = "in", 
           plot = pl, 
           width = 5, 
           height = 4)
  }
  
  return(ktest_all )
}

test_distributions_per_class <- function(data_adj, 
                                         data_or, 
                                         outfold = NULL, 
                                         save_plot = FALSE, 
                                         show_plot = TRUE, 
                                         save_ext = "png"){
  
  CL_names <- colnames(data_adj)[-(1:13)]
  # remove AnchorCombinations
  data_adj <- data_adj %>%
    dplyr::filter(!Note1 %in% "AnchorCombinations")
  
  data_or <- data_or %>%
    dplyr::filter(!Note1 %in% "AnchorCombinations")
  
  class_names <- unique(data_adj$Note1)
  ktest_adj <- ktest_or <- list()
  for (i in 1:length(CL_names)) {
    CL <- CL_names[i]
    ktest_adj[[i]] <- sapply(class_names, function(x) kruskal.test(x = data_adj[data_adj$Note1 %in% x, CL], g = data_adj[data_adj$Note1 %in% x, "lib"])$p.value)
    ktest_adj[[i]] <- data.frame(pvalue = unname(ktest_adj[[i]]), class = names(ktest_adj[[i]]), CL = CL)
    ktest_or[[i]] <- sapply(class_names, function(x) kruskal.test(x = data_or[data_or$Note1 %in% x, CL], g = data_or[data_or$Note1 %in% x, "lib"])$p.value)
    ktest_or[[i]] <- data.frame(pvalue = unname(ktest_or[[i]]), class = names(ktest_or[[i]]), CL = CL)
  }
  ktest_adj <- do.call(rbind, ktest_adj)
  ktest_or <- do.call(rbind, ktest_or)
  ktest_all <- full_join( ktest_adj, ktest_or, 
                          by = c("CL", "class"), 
                          suffix = c("_adjusted", "_original")) %>%
    dplyr::mutate(log10_pvalue_adjusted = -log10(pvalue_adjusted + .Machine$double.xmin), 
                  log10_pvalue_original = -log10(pvalue_original + .Machine$double.xmin))
  
  pl <- ggplot(ktest_all, aes(x = log10_pvalue_adjusted,
                              y = log10_pvalue_original,
                              col = CL,
                              label = CL)) +
    facet_wrap(.~class, ncol = 2) +
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", col = "red") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "red") +
    # geom_text_repel(size = 3, max.overlaps = Inf) +
    theme_bw() + 
    theme(legend.position = "right", legend.title = element_blank()) + 
    xlab("Corrected") + 
    ylab("Raw") + 
    ggtitle("-log10 P-value from Kruskal-Wallis test")
  
  if (show_plot) {
    print(pl)
  }
  
  if (save_plot) {
    ggsave(filename = sprintf("%skruskaltest_perclass_perCL.%s", outfold, save_ext), 
           units = "in", 
           plot = pl, 
           width = 6, 
           height = 8)
  }
  
  return(ktest_all)
  
}

plot_library_genes <- function(data_adj, 
                               data_or, 
                               essential_genes, 
                               outfold = NULL, 
                               show_plot = TRUE,
                               save_plot = FALSE, 
                               save_ext = "png"){
  
  data_annot_lib <- data_adj %>%
    dplyr::filter(Note1 %in% c("LibrarySingletons", "LibraryCombinations")) %>%
    dplyr::mutate(Gene1_essential = Gene1 %in% essential_genes, 
                  Gene1_or_Gene2_essential = Gene1 %in% essential_genes | Gene2 %in% essential_genes) 
  
  data_adj_lib <- data_adj %>%
    dplyr::filter(Note1 %in% c("LibrarySingletons", "LibraryCombinations"))
  data_or_lib <- data_or %>%
    dplyr::filter(Note1 %in% c("LibrarySingletons", "LibraryCombinations"))
  CL_names <- colnames(data_adj_lib)[-(1:13)]
  other_names <- colnames(data_adj_lib)[1:13]
  
  tmp <- lapply(CL_names, function(x) data_adj_lib[, c(other_names, x)] %>% 
                  dplyr::rename(logFC = !!(x)) %>% 
                  dplyr::mutate(CL = x))
  data_adj_lib <- do.call(rbind, tmp)
  
  tmp <- lapply(CL_names, function(x) data_or_lib[, c(other_names, x)] %>% 
                  dplyr::rename(logFC = !!(x)) %>% 
                  dplyr::mutate(CL = x))
  data_or_lib <- do.call(rbind, tmp)
  
  pl1 <- ggplot(data_adj_lib, 
                aes(x = lib, y = logFC, fill = lib)) + 
    geom_boxplot(outlier.size = 1) + 
    theme_bw() + 
    xlab("") +
    ggtitle(sprintf("%s: Corrected", "All CLs")) + 
    theme(legend.position = "none")
  
  pl2 <- ggplot(data_or_lib, 
                aes(x = lib, y = logFC, fill = lib)) + 
    geom_boxplot(outlier.size = 1) + 
    theme_bw() + 
    xlab("") +
    ggtitle(sprintf("%s: Raw", "All CLs")) + 
    theme(legend.position = "none")
  
  pl_lib <- ggpubr::ggarrange(plotlist = list(pl2, pl1),
                              ncol = 2)
  
  data_lib_perc <- data_annot_lib %>%
    dplyr::group_by(lib) %>%
    dplyr::summarise(frac_Gene1_ess = sum(Gene1_essential)/length(Gene1_essential), 
                     n_Gene1_ess = sum(Gene1_essential),
                     frac_Gene1_or_Gene2_ess = sum(Gene1_or_Gene2_essential)/length(Gene1_or_Gene2_essential), 
                     n_Gene1_or_Gene2_ess = sum(Gene1_or_Gene2_essential), 
                     n_tot = length(Gene1_essential))
  
  pl_frac <- ggplot(data_lib_perc , aes(x = frac_Gene1_ess,
                                        y = n_Gene1_ess,
                                        col = lib,
                                        label = lib)) +
    geom_point(size = 2) +
    # geom_line(alpha = 0.5) + 
    geom_text_repel(size = 3, max.overlaps = Inf) +
    guides(color = "none", fill = "none") +
    theme_bw() + 
    xlab("Fraction of pairs with Gene1 essential") + 
    ylab("N. of pairs with Gene1 essential")
  
  if (show_plot) {
    print(pl_lib)
    print(pl_frac)
  }
  
  if (save_plot) {
    ggsave(filename = sprintf("%sALLCLs_distr_perLib_LibraryGenes.%s", outfold, save_ext), 
           units = "in", 
           plot = pl_lib, 
           width = 5, 
           height = 5)
    
    ggsave(filename = sprintf("%sFraction_Gene1_ess_perLib_LibraryGenes.%s", outfold, save_ext), 
           units = "in", 
           plot = pl_frac, 
           width = 3, 
           height = 3)
  }
  
  return(list(logFC_or = data_or_lib, 
              logFC_adj = data_adj_lib, 
              perc = data_lib_perc))
  
}






