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
    ggsave(filename = sprintf("%sComBat_param_dist_libraries.png", outfold), 
           units = "in", 
           plot = pl_lib, 
           width = 4, 
           height = 4)
    
    ggsave(filename = sprintf("%sComBat_param_dist_GPclass.png", outfold), 
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


