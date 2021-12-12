basedir <- "path/do/basedir"

suppressPackageStartupMessages({
  library("TCA")
  library("EpiDISH")
  library("matrixStats")
  require("GEOquery")
  library("R.utils")
  require("data.table")
  library("pracma")
  library("stats")
  library("stringr")
  library("ggplot2")
  library("ggpubr")
  library("grid")
  library("hrbrthemes")
  library("extrafont")
  extrafont::loadfonts()
})


# A function for downloading the Liu et al. and Hannum et al. Data
# data_dir - destination directory for RData files with the data
# hannum_smk_status_path - Path to hannum_smoking_status.txt that must be downloaded separately.
# Returns a list of filenames for these RData files
prep_data <- function(data_dir, hannum_smk_status_path){
  options(timeout=10000)
  file_name1 <- paste(data_dir,"liu.RData",sep="/")
  file_name1_processed <- paste(data_dir,"liu.processed.RData",sep="/")
  if (!file.exists(file_name1)){
    # Download the Liu et al. data
    # note that the GEO series matrix of the Liu et al. data was updates on 2021-09-20; the new version now includes only ~140k CpGs. We therefore download the full processed version of the data.
    gse <- GEOquery::getGEO("GSE42861", destdir = data_dir, GSEMatrix = TRUE)
    download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/suppl/GSE42861_processed_methylation_matrix.txt.gz", file.path(data_dir,"GSE42861_processed_methylation_matrix.txt.gz"))
    gunzip(file.path(data_dir,"GSE42861_processed_methylation_matrix.txt.gz"))
    X <- fread(file.path(data_dir,"GSE42861_processed_methylation_matrix.txt"))
    cpg_names <- X$ID_REF
    sample_names <- colnames(X)[2:(ncol(X)-1)] # the last column is pvals
    X <- as.matrix(X[,2:(ncol(X)-1)])
    colnames(X) <- sample_names
    rownames(X) <- cpg_names
    str_split(Biobase::pData(gse[[1]])[["supplementary_file"]],"_")
    ids_map <- rownames(Biobase::pData(gse[[1]]))
    l <- str_split(Biobase::pData(gse[[1]])[["supplementary_file"]],"_")
    ids_map_names <- character(length(l))
    for (i in 1:length(l)){
      ids_map_names[i] <- paste(l[[i]][2],l[[i]][3],sep="_")
    }
    names(ids_map) <- ids_map_names
    X <- X[,names(ids_map)]
    colnames(X) <- as.character(ids_map)
    
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
    rownames(cov.liu) <- rownames(Biobase::pData(gse[[1]]))[keep]
    rownames(smk.liu) <- rownames(cov.liu)
    # methylaion data
    X.liu <- X[,rownames(smk.liu)]
    
    # estimate cell-type proportions using a reference-based approach 
    # methylation reference
    ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
    W.liu <- EpiDISH::epidish(X.liu, ref)$estF
    
    liu <- list(X=X.liu,
                smk=smk.liu,
                cov=cov.liu,
                W=W.liu)
    gse.file <- paste0(data_dir, "/GSE42861_series_matrix.txt.gz")
    other.gse.file <- paste0(data_dir, "/GPL13534.soft")
    file.remove(gse.file)
    file.remove(other.gse.file)
    save(liu, file=file_name1)
    rm(gse, X.liu, smk.liu, cov.liu, W.liu, liu)
    
  }
  
  file_name2 <- paste(data_dir,"hannum.RData",sep="/")
  file_name2_processed <- paste(data_dir,"hannum.processed.RData",sep="/")
  if (!file.exists(file_name2)){
    # Download the Hannum et al. data 
    gse <- GEOquery::getGEO("GSE40279", destdir = data_dir, GSEMatrix = TRUE)
    # smoking status; consider never-smokers [0], ex-smokers [1] and current smokers [2] (refer to occasional smokers as smokers)
    smk.hannum <- read.table(hannum_smk_status_path, header = T, row.names = 1, sep=",")
    # remove samples with NA values for smoking
    keep <- setdiff(1:nrow(smk.hannum), which(is.na(smk.hannum)))
    keep.ids <- rownames(smk.hannum)[keep]
    smk.hannum <- smk.hannum[keep.ids,,drop=F]
    smk.hannum[smk.hannum=="Never"] <- 0
    smk.hannum[smk.hannum=="Past Smoker"] <- 1
    smk.hannum[smk.hannum=="Current Smoker"] <- 2
    smk.hannum <- as.numeric(smk.hannum$smk)
    names(smk.hannum) <- keep.ids
    # covariates
    cov.hannum <- cbind(as.numeric(as.factor(Biobase::pData(gse[[1]])[["age (y):ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["gender:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["plate:ch1"]])), as.numeric(as.factor(Biobase::pData(gse[[1]])[["ethnicity:ch1"]])))
    rownames(cov.hannum) <- Biobase::pData(gse[[1]])[["geo_accession"]]
    colnames(cov.hannum) <- c("age", "gender", "plate","ethnicity")
    cov.hannum <- cov.hannum[keep.ids,]
    # methylaion data
    X.hannum <- Biobase::exprs(gse[[1]])[,keep.ids]
    # remove sites with missing values 
    X.hannum <- X.hannum[rowSums(is.na(X.hannum)) == 0,]
    
    names(smk.hannum) <- colnames(X.hannum)
    rownames(cov.hannum) <- colnames(X.hannum)
    
    # estimate cell-type proportions using a reference-based approach 
    ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
    W.hannum <- EpiDISH::epidish(X.hannum, ref)$estF
    
    hannum <- list(X=X.hannum,
                   smk=smk.hannum,
                   cov=cov.hannum,
                   W=W.hannum)
    gse.file <- paste0(data_dir, "/GSE40279_series_matrix.txt.gz")
    file.remove(gse.file)
    file.remove(other.gse.file)
    save(hannum, file=file_name2)
    rm(gse, X.hannum, smk.hannum, cov.hannum, W.hannum)
    
  }
  if (!file.exists(file_name1_processed) | !file.exists(file_name2_processed)){
    
    load(file_name1)
    load(file_name2)
    
    # calculate PCs from low variance probes , to be treated as control probes (as in Lenhe et al. 2015); here, since we don't work with IDAT files and therefore don't have actual control probes, we use sites with the least variation in the data as control probes (as in Rahmani et al. 2019).
    p <- 10000
    site.variances1 <- matrixStats::rowVars(liu$X)
    names(site.variances1) <- rownames(liu$X)
    low.var.sites1 <- names(head(sort(site.variances1), p))
    low.var.pca1 <- prcomp(t(liu$X[low.var.sites1,]), center=TRUE, scale=TRUE, rank=20)
    liu$ctrl_pcs <- low.var.pca1$x
    
    site.variances2 <- matrixStats::rowVars(hannum$X)
    names(site.variances2) <- rownames(hannum$X)
    low.var.sites2 <- names(head(sort(site.variances2), p))
    low.var.pca2 <- prcomp(t(hannum$X[low.var.sites2,]), center=TRUE, scale=TRUE, rank=20)
    hannum$ctrl_pcs <- low.var.pca2$x
    
    # remove polymorphic or cross-reactive sites and non-autosomal sites
    nonspecific_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/nonspecific_probes.txt")[,1]
    XY_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/HumanMethylationSites_X_Y.txt")[,1]
    polymorphic_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/polymorphic_cpgs.txt")[,1]
    
    # remove sites with very low variance that are unlikely to exhibit biological signal
    low_var_sites1 <- names(which(site.variances1<0.001))
    low_var_sites2 <- names(which(site.variances2<0.001))
    
    # in the set of cpgs to exclude do not include known smoking cpgs
    exclude1 <- union(nonspecific_probes,union(XY_probes,union(low_var_sites1,polymorphic_probes)))
    
    exclude2 <- union(nonspecific_probes,union(XY_probes,union(low_var_sites2,polymorphic_probes)))
    
    keep <- intersect(setdiff(rownames(liu$X),exclude1), setdiff(rownames(hannum$X),exclude2))
    
    liu$X <- liu$X[keep,]
    hannum$X <- hannum$X[keep,]
    
    liu$cov <- cbind(liu$cov,liu$smk)
    liu$smk <- NULL
    W_colnames <- c("Gran",setdiff(colnames(liu$W), c("Neutro","Eosino")))
    liu$W <- cbind(rowSums(liu$W[,c("Neutro","Eosino")]), liu$W[,setdiff(colnames(liu$W), c("Neutro","Eosino"))])
    colnames(liu$W) <- W_colnames
    
    cov_colnames <- c(colnames(hannum$cov),"smoking") 
    hannum$cov <- cbind(hannum$cov,hannum$smk)
    colnames(hannum$cov) <- cov_colnames
    hannum$smk <- NULL
    W_colnames <- c("Gran",setdiff(colnames(hannum$W), c("Neutro","Eosino")))
    hannum$W <- cbind(rowSums(hannum$W[,c("Neutro","Eosino")]), hannum$W[,setdiff(colnames(hannum$W), c("Neutro","Eosino"))])
    colnames(hannum$W) <- W_colnames
    
    save(liu, file=file_name1_processed)
    save(hannum, file=file_name2_processed)
    
  }
  return(list(file_name1_processed, file_name2_processed, file_name1, file_name2))
}


plot_consistency_analsyis <- function(results1, results2, plot_title, nums = 1:100, interval_for_plot = seq(10,100,10)){
  
  score.tca <- numeric(length(nums))
  score.tcareg <- numeric(length(nums))
  score.celldmc <- numeric(length(nums))
  results1$tca.order <- order(results1$tca)
  results2$tca.order <- order(results2$tca)
  results1$tcareg.order <- order(results1$tcareg)
  results2$tcareg.order <- order(results2$tcareg)
  results1$celldmc.order <- order(results1$celldmc)
  results2$celldmc.order <- order(results2$celldmc)
  for (i in 1:length(nums)){
    score.tca[i] <- length(intersect(results1$tca.order[1:nums[i]],results2$tca.order[1:nums[i]]))/nums[i]
    score.tcareg[i] <- length(intersect(results1$tcareg.order[1:nums[i]],results2$tcareg.order[1:nums[i]]))/nums[i]
    score.celldmc[i] <- length(intersect(results1$celldmc.order[1:nums[i]],results2$celldmc.order[1:nums[i]]))/nums[i]
  }
  
  df <- data.frame(x = interval_for_plot,
                   score = c(score.tca[interval_for_plot],
                             score.celldmc[interval_for_plot],
                             tcareg = score.tcareg[interval_for_plot]),
                   Method = c(rep("TCA X|Y",length(interval_for_plot)),
                               rep("CellDMC X|Y",length(interval_for_plot)),
                              rep("TCA Y|X",length(interval_for_plot))))
  
  res <- list()              
  res$plot <- ggplot(data = df, aes(x = x, y = score, group = Method)) +
    geom_line(aes(linetype=Method, color=Method), size=0.75) +
    theme_ipsum(axis_title_just="c") +
    xlab("# Top significant CpGs") + ylab("Validation rate") +
    ggtitle(plot_title) +
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
    theme(legend.position="bottom") +
    theme(plot.title = element_text(hjust = 0.5),legend.margin=margin(t = 0, unit='cm'),axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), legend.title=element_blank())
  
  return(res)
}

extract_pvals <- function(liu_results, hannum_results, cell_types){
  liu_pvals <- list()
  hannum_pvals <- list()
  liu_pvals$celldmc <- cbind(liu_results$celldmc$coe$Gran$p,liu_results$celldmc$coe$CD4T$p,
                             liu_results$celldmc$coe$CD8T$p,liu_results$celldmc$coe$Mono$p,
                             liu_results$celldmc$coe$B$p,liu_results$celldmc$coe$NK$p)
  hannum_pvals$celldmc <- cbind(hannum_results$celldmc$coe$Gran$p,hannum_results$celldmc$coe$CD4T$p,
                                hannum_results$celldmc$coe$CD8T$p,hannum_results$celldmc$coe$Mono$p,
                                hannum_results$celldmc$coe$B$p,hannum_results$celldmc$coe$NK$p)
  liu_pvals$tca <- liu_results$tca$gammas_hat_pvals[,paste(cell_types,"y",sep=".")]
  hannum_pvals$tca <- hannum_results$tca$gammas_hat_pvals[,paste(cell_types,"y",sep=".")]
  liu_pvals$tcareg <- liu_results$tcareg$pvals
  hannum_pvals$tcareg <- hannum_results$tcareg$pvals  
  return(list(liu_pvals, hannum_pvals))
}


run_ewas <- function(X, W, y, C1, C2){
  results <- list()
  print("celldmc...")
  results$celldmc <- CellDMC(X, y, W, cov.mod = cbind(C1,C2))
  print("tca...")
  results$tca <- tca(X = X, W = W, C1 = cbind(C1, y), C2 = C2, verbose = TRUE)
  print("tcareg...")
  tca.mdl <- tca(X,W,C1=C1,C2=C2, verbose = TRUE)
  results$tcareg <- tcareg(X, tca.mdl, y ,C3 = cbind(C1,C2), test = "marginal_conditional",save_results = FALSE, log_file = NULL, verbose = TRUE)
  return(results)
}



data_dir <- file.path(basedir,"data/") 
results_dir <- file.path(basedir,"results/")
hannum_smk_status_path <- file.path(basedir,"assets/hannum_smoking_status.txt") # hannum_smoking_status.txt needs to be obtained from the authors of Hannum et al.
outfile <- file.path(results_dir,"ewas_results.pdf")

dir.create(data_dir, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)

# download and prepare data for the analysis
data.paths <- prep_data(data_dir, hannum_smk_status_path)

liu_results_file_age <- paste(results_dir,"liu_age_results.RData",sep="/")
liu_results_file_sex <- paste(results_dir,"liu_sex_results.RData",sep="/")
hannum_results_file_age <- paste(results_dir,"hannum_age_results.RData",sep="/")
hannum_results_file_sex <- paste(results_dir,"hannum_sex_results.RData",sep="/")

load(data.paths[[1]]) # Load Liu
load(data.paths[[2]]) # Load Hannum
cell_types <- paste(colnames(liu$W))

# run the analysis
if (!file.exists(liu_results_file_age)){
  liu_results <- run_ewas(X = liu$X, W=liu$W, y = liu$cov[,"age"], C1 = liu$cov[, c("disease","gender","smoking")], C2 = liu$ctrl_pcs)
  save(liu_results, file=liu_results_file_age)
}
if (!file.exists(hannum_results_file_age)){
  hannum_results <- run_ewas(X = hannum$X, W=hannum$W, y = hannum$cov[,"age"], C1 = hannum$cov[, c("ethnicity","gender","smoking")], C2 = hannum$ctrl_pcs)
  save(hannum_results, file=hannum_results_file_age)
}

if (!file.exists(liu_results_file_sex)){
  liu_results <- run_ewas(X = liu$X, W=liu$W, y = liu$cov[,"gender"], C1 = liu$cov[, c("disease","age","smoking")], C2 = liu$ctrl_pcs)
  save(liu_results, file=liu_results_file_sex)
}
if (!file.exists(hannum_results_file_sex)){
  hannum_results <- run_ewas(X = hannum$X, W=hannum$W, y = hannum$cov[,"gender"], C1 = hannum$cov[, c("ethnicity","age","smoking")], C2 = hannum$ctrl_pcs)
  save(hannum_results, file=hannum_results_file_sex)
}


# summarize age analysis
load(liu_results_file_age)
load(hannum_results_file_age)
pvals <- extract_pvals(liu_results, hannum_results, cell_types)
res_age <- plot_consistency_analsyis(pvals[[1]], pvals[[2]], plot_title = "Age")

# summarize sex analysis
load(liu_results_file_sex)
load(hannum_results_file_sex)
pvals <- extract_pvals(liu_results, hannum_results, cell_types)
res_sex <- plot_consistency_analsyis(pvals[[1]], pvals[[2]], plot_title = "Sex")

# plot results
plt <- ggarrange(res_age$plot, res_sex$plot, 
                     ncol = 2, nrow = 1, labels = c("a","b"))
ggsave(outfile, plot = plt, width = 6.5, height = 3.5)
