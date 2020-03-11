process_meth_data <- function(X){
  # remove sites with missing values 
  X <- X[rowSums(is.na(X)) == 0,]
  # remove consistently hyper or hypo methylated cpgs
  means <- rowMeans(X)
  X <- X[(means >= 0.1) & (means <= 0.9), ]
}

#' Download and organize Liu et al. and Hannum et al. Data
#'
#' @param data_path Directory for RData files that will store data
#' @param hannum_smk_status_path Path to hannum_smoking_status.txt
#'        that must be downloaded separately.
#' @return List of filenames for these RData files
prep_smoking_data <- function(data_path, hannum_smk_status_path){
  file_name1 <- paste(data_path,"liu.RData",sep="/")
  if (!file.exists(file_name1)){
    # Download the Liu et al. data 
    gse <- GEOquery::getGEO("GSE42861", destdir = data_path, GSEMatrix = TRUE)
    # smoking status; consider never-smokers [0], ex-smokers [1] and current smokers [2] (refer to occasional smokers as smokers)
    smk.liu <- Biobase::pData(gse[[1]])[["smoking status:ch1"]]
    # two samples have NA values for smoking; remove them from the analysis
    keep <- setdiff(1:length(smk.liu), which(smk.liu == "na"))
    smk.liu <- smk.liu[keep]
    smk.liu[which(smk.liu == "occasional")] <- "current"
    smk.liu <- 3-as.matrix(as.numeric(as.factor(smk.liu)))
    colnames(smk.liu) <- "smoking"
    # covariates
    cov.liu <- cbind(as.numeric(as.factor(Biobase::pData(gse[[1]])[["disease state:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["gender:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["age:ch1"]])))
    colnames(cov.liu) <- c("disease", "gender", "age")
    cov.liu <- cov.liu[keep,]
    # methylaion data
    X.liu <- process_meth_data(Biobase::exprs(gse[[1]])[,keep])
    rownames(smk.liu) <- colnames(X.liu)
    rownames(cov.liu) <- colnames(X.liu)
    
    # estimate cell-type proportions usign a reference-based approach 
    # methylation reference
    ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
    W.liu <- EpiDISH::epidish(X.liu, ref)$estF
    
    #save(X.liu, smk.liu, cov.liu, W.liu, ref.liu, file = file_name1)
    liu <- list(X=X.liu,
                smk=smk.liu,
                cov=cov.liu,
                W=W.liu)
    save(liu, file=file_name1)
    rm(gse, X.liu, smk.liu, cov.liu, ref.liu)
  }
  
  file_name2 <- paste(data_path,"hannum.RData",sep="/")
  if (!file.exists(file_name2)){
    # Download the Hannum et al. data 
    gse <- GEOquery::getGEO("GSE40279", destdir = data_path, GSEMatrix = TRUE)
    
    ## TODO download from github instead
    # smoking status; consider never-smokers [0], ex-smokers [1] and current smokers [2] (refer to occasional smokers as smokers)
    smk.hannum <- read.table(hannum_smk_status_path, header = T, row.names = 1, sep=",")
    # remove samples with NA values for smoking
    keep <- setdiff(1:nrow(smk.hannum), which(is.na(smk.hannum)))
    keep.ids <- rownames(smk.hannum)[keep]
    smk.hannum <- smk.hannum[keep.ids,,drop=F]
    smk.hannum$smk <- as.numeric(smk.hannum$smk)
    smk.hannum$smk[which(smk.hannum$smk == 2)] <- 0
    smk.hannum$smk[which(smk.hannum$smk == 1)] <- 2
    smk.hannum$smk[which(smk.hannum$smk == 3)] <- 1
    colnames(smk.hannum) <- "smoking"
    # covariates
    cov.hannum <- cbind(as.numeric(as.factor(Biobase::pData(gse[[1]])[["age (y):ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["gender:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["plate:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["ethnicity:ch1"]])))
    rownames(cov.hannum) <- Biobase::pData(gse[[1]])[["geo_accession"]]
    colnames(cov.hannum) <- c("age", "gender", "plate","ethnicity")
    cov.hannum <- cov.hannum[keep.ids,]
    # methylaion data
    X.hannum <- process_meth_data(Biobase::exprs(gse[[1]])[,keep.ids])
    rownames(smk.hannum) <- colnames(X.hannum)
    rownames(cov.hannum) <- colnames(X.hannum)
    
    # estimate cell-type proportions usign a reference-based approach 
    ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
    W.hannum <- EpiDISH::epidish(X.hannum, ref)$estF
    
    hannum <- list(X=X.hannum,
                   smk=smk.hannum,
                   cov=cov.hannum,
                   W=W.hannum)
    save(hannum, file=file_name2)
    #save(X.hannum, smk.hannum, cov.hannum, W.hannum, ref.hannum, file = file_name2)
    rm(gse, X.hannum, smk.hannum, cov.hannum, ref.hannum)
    
  }
  
  return(list(file_name1, file_name2))
  
}

#' Smoking analysis of methylation data
#'
#' Performs cell-type-specific analysis of blood methylation data in the context of smoking.
#' Provides p-values estimated by CellDMC and TCA while correcting for given covariates
#' as well as ReFACTor components. Combines cell types into myeloid and lymphoid compartments
#' for analysis.
#'
#' @param X Methylation values in matrix with sites as rows and samples as columns
#' @param smk Smoking status encoded as 0 (non-smoker), 1 (ex-smoker), or 2 (smoker) in
#'            a matrix with rownames as samples and 1 column entitled \code{smoking}
#' @param covs Covariates in matrix format. Rownames are samples and column names are
#'                   covariate names
#' @param W Cell fractions in a matrix with cell types as column names and samples as
#'          row names. Expects cell types: "Neutro","Eosino","Mono","CD4T","CD8T","B","NK"
#' @param smoking_cpgs List of CpG names to analyze for associations
#' @param dataset.name Name of dataset for naming output variables.
#' @return A list of results. Each slot contains p-values output by each method
run_smoking_analysis <- function(X, smk, covs, W, smoking_cpgs, dataset.name){
  # combine cell types into myeloid and lymphoid compartments
  W.2 <- cbind(rowSums(W[,c("Neutro","Eosino","Mono")]),
               rowSums(W[,c("CD4T","CD8T","B","NK")]))
  colnames(W.2) <- c("MYE","LYM")
  
  # Compute the ReFACTor components of the data as additional covariates
  refactor.components <- TCA::refactor(X, k=6, rand_svd=T, C=covs)
  
  # run celldmc and tca and extract p-values
  celldmc.res <- CellDMC(X[smoking_cpgs,], smk, W.2,
                         cov.mod = cbind(covs,refactor.components$scores))
  celldmc.pvals <- list(celldmc.res$coe$MYE$p, celldmc.res$coe$LYM$p)
  
  # run tca and extract p-values
  tca.mdl <- tca(X = X[smoking_cpgs,], W = W.2,
                 C1 = cbind(covs, smk),
                 C2 = refactor.components$scores)
  tca.pvals <- list(tca.mdl$gammas_hat_pvals[,"MYE.smoking"], 
                    tca.mdl$gammas_hat_pvals[,"LYM.smoking"])
  tca.pvals.joint <- list(tca.mdl$gammas_hat_pvals.joint[,"smoking"])
  
  res.list <- list()
  res.list[[sprintf("celldmc.pvals.%s", dataset.name)]] <- celldmc.pvals
  res.list[[sprintf("tca.pvals.%s", dataset.name)]] <- tca.pvals
  res.list[[sprintf("tca.pvals.joint.%s", dataset.name)]] <- tca.pvals.joint
  return(res.list)
}



plot_smoking_analysis <- function(smoking.res, outfile){
  p.celldmc <- plot_smoking_heatmap(smoking.res$smoking_cpgs, smoking.res$celldmc.pvals.liu, smoking.res$celldmc.pvals.hannum, smoking.res$cell_types, "CellDMC")
  p.tca <- plot_smoking_heatmap(smoking.res$smoking_cpgs, smoking.res$tca.pvals.liu, smoking.res$tca.pvals.hannum, smoking.res$cell_types, "TCA")
  p.tca.joint <- plot_smoking_heatmap(smoking.res$smoking_cpgs, smoking.res$tca.pvals.joint.liu, smoking.res$tca.pvals.joint.hannum, c("joint"), "TCA - joint test")
  p.final <- ggarrange(p.celldmc, p.tca, p.tca.joint, ncol = 3, nrow = 1, labels = c("a","b","c"))
  ggsave(outfile, plot = p.final, width = 8.5, height = 3)
  return(p.final)
}


plot_smoking_heatmap <- function(cpgs, pvals.liu, pvals.hannum, cell_types, method_title){
  data <- data.frame(matrix(ncol = 4, nrow = 2*length(cell_types)*length(cpgs)))
  colnames(data) <- c("cpg","cell_type","neg_log_pval","dataset")
  for (i in 1:length(cell_types)){
    data[((i-1)*length(cpgs)+1):(i*length(cpgs)),] <- cbind(cpgs,
                                                            rep(c(cell_types[[i]])),
                                                            -log10(pvals.hannum[[i]]),
                                                            rep(c("Hannum et al."),
                                                                length(cpgs)))
  }
  for (i in 1:length(cell_types)){
    data[(length(cell_types)*length(cpgs)+(i-1)*length(cpgs)+1):(length(cell_types)*length(cpgs)+i*length(cpgs)),] <- cbind(cpgs,
                                                                                                                            rep(c(cell_types[[i]])),
                                                                                                                            -log10(pvals.liu[[i]]),
                                                                                                                            rep(c("Liu et al."),
                                                                                                                                length(cpgs)))
  }
  data$neg_log_pval <- as.numeric(data$neg_log_pval)
  data$cpg <- factor(data[,1],rev(data[1:length(cpgs),1]))
  axis.text.x <- if (length(cell_types)-1) element_text(colour = c("#009E73","#D55E00" )) else axis.text.x = NULL
  plt <- ggplot(data, aes(cell_type, cpg, fill= neg_log_pval)) +  geom_tile() + 
    scale_fill_gradient(low = "#808080", high = "#FFFF00", na.value = "#FFFF00", limits = c(0,15)) +
    geom_text(aes(label = round(data$neg_log_pval,1))) + 
    geom_hline(aes(yintercept=2.5), size = 1) +
    facet_wrap(~dataset) +
    ggtitle(method_title) +
    ylab(NULL) +
    xlab(NULL) +
    labs(fill="-log10(P)") +
    theme(axis.text.y = element_text(colour = c("#009E73","#009E73","#D55E00","#D55E00",
                                                "#D55E00","#D55E00","#D55E00")),
          axis.text.x = axis.text.x, plot.title = element_text(hjust = 0.5),
          strip.background = element_rect(fill="white"),
          panel.background = element_rect(fill = "white", colour = "white"),
          legend.position="bottom", strip.text = element_text(size = 10,face = "bold")) 
  return(plt)
}

#' Replicate smoking analysis of Liu et al. and Hannum et al. datasets
#'
#' @param liu List generated by get_liu_data() or provided in CellTypeSpecificMethylationData
#'            package under the same name
#' @param hannum List generated by get_hannum_data() or provided in CellTypeSpecificMethylationData
#'               package under the same name
#' @param results_dir Directory to save plot of results
#' @param plot_type Extension for saving plot graphics, such as pdf or png
#' @return A list of results. Slot \code{results} contains a list of p-values generated by each method.
#'         Slot \code{plot} contains plot of results.
#' @export
smoking_analysis <- function(liu, hannum, results_dir, plot_type){
  # Smoking cpgs from Su et al. that were considered by Jing et al.
  smoking_cpgs <- c("cg05575921","cg21566642","cg09935388","cg06126421",
                    "cg03636183","cg19859270","cg09099830")
  
  # Cell types will be combined into myeloid and lymphoid compartments
  cell.types <- c("MYE","LYM")
  
  # Liu et al. analysis
  liu.smoking.res <- run_smoking_analysis(X = liu$X,
                                          smk = liu$smk, 
                                          covs = liu$cov, 
                                          W = liu$W, 
                                          smoking_cpgs = smoking_cpgs, 
                                          dataset.name = "liu")
  # Hannum et al. analysis
  hannum.smoking.res <- run_smoking_analysis(X = hannum$X,
                                             smk = hannum$smk, 
                                             covs = hannum$cov, 
                                             W = hannum$W, 
                                             smoking_cpgs = smoking_cpgs, 
                                             dataset.name = "hannum")
  
  smoking.res <- c(liu.smoking.res,
                   hannum.smoking.res,
                   list("cell_types" = cell.types,
                        "smoking_cpgs" = smoking_cpgs))
  outfile <- paste(results_dir, "/smoking_analysis.", plot_type,
                   sep="")
  smoking.plot <- plot_smoking_analysis(smoking.res, outfile)
  return(list(smoking.res=smoking.res,
              smoking.plot=smoking.plot))
}
