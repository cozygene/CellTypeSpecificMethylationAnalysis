library(EpiDISH)
library(TCA)
library(ggplot2)
library(plyr)
library(cowplot)

#' Compares cell fraction estimates from EpiDISH and TCA under several models
#'
#' @param X Matrix of methylation values with column names as samples and row names as sites
#' @param W.facs Matrix of FACS derived cell type fractions with column names as cell types
#'               and row names as samples
#' @param C1 Matrix of source-specific covariates for TCA with column names as features and
#'           row names as samples
#' @param C2 Matrix of mixture-level covariates for TCA with column names as features and
#'           row names as samples
#' @param random_seed Random seed for replication
#' @param verbose Boolean for printing out progress
#' @return A list of results. Slot performance.metrics contains a dataframe where each
#'         row indicates either the PCC or RMSE of a method for estimating a specific
#'         cell type fraction. Slot prop.estimates contains a list of fractions estimated
#'         by each method.
compare_cell_fraction_estimates <- function(X, W.facs, C1=NULL, C2=NULL, random_seed=1000, verbose=TRUE){
  set.seed(random_seed)
  cell.types <- c('B','NK','CD4T','CD8T','Mon', 'Neu')
  if (all.equal(sort(colnames(W.facs)),sort(cell.types)) != TRUE){
    stop("W.facs should have following column names: ",
         "\'B\',\'NK\',\'CD4T\',\'CD8T\',\'Mon\', \'Neu\'")
  }
  samples <- colnames(X)
  ref.sites <- rownames(EpiDISH::centDHSbloodDMC.m)
  available.ref.sites <- intersect(ref.sites,rownames(X))
  if (verbose){
    message(sprintf("%i reference sites present.", length(available.ref.sites)))
  }
  cell.types <- c('B','NK','CD4T','CD8T','Mono','Neutro')
  
  # EpiDISH estimates
  epidish.est.props <- EpiDISH::epidish(beta.m = X,
                                        ref.m=centDHSbloodDMC.m[available.ref.sites,cell.types],
                                        method = "RPC")$estF
  # Match cell type names
  colnames(epidish.est.props) <- plyr::mapvalues(colnames(epidish.est.props),
                                                 from = c('B','NK','CD4T','CD8T','Mono','Neutro'),
                                                 to = c('B','NK','CD4T','CD8T','Mon', 'Neu'))
  cell.types <- c('B','NK','CD4T','CD8T','Mon', 'Neu')
  epidish.est.props <- epidish.est.props[,cell.types]
  
  # Noisy EpiDISH estimates
  min.w <- 0.01
  noise <- sapply(cell.types, function(cell.type){
    noise.sd <- mean(epidish.est.props[,cell.type])**2
    noise.dir <- sample(c(-1,1), 1, prob=c(0.5, 0.5))
    return(abs(rnorm(nrow(epidish.est.props), sd = noise.sd)) * noise.dir)
  })
  noisy.epidish.props <- epidish.est.props + noise
  noisy.epidish.props[noisy.epidish.props <= 0] <- min.w
  noisy.epidish.props <- noisy.epidish.props/rowSums(noisy.epidish.props)
  
  # Get low variance site PCs
  message("Starting low variance site PCA")
  site.variances <- matrixStats::rowVars(X)
  names(site.variances) <- rownames(X)
  p <- 1000
  low.var.sites <- names(head(sort(site.variances), p))
  low.var.pca <- prcomp(t(X[low.var.sites,]), center=TRUE, scale=TRUE, rank=10)
  # Covariate data for refactor
  cov.data <- low.var.pca$x
  if (!is.null(C2)){
    cov.data <- cbind(cov.data, C2)
  }
  
  # Get refactor sites
  # run refactor as follows, according to the guidlines in https://glint-epigenetics.readthedocs.io/en/latest/tissueheterogeneity.html#refactor
  # remove polymorphic or cross-reactive sites and non-autosomal sites
  # we also consider control-prob based PCs as covaraites (as in Lenhe et al. 2015); here, since we don't work with IDAT files and therefore don't have actual control probes, we use sites with the least variation in the data as control probes (as in Rahmani et al. 2019).
  nonspecific_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/nonspecific_probes.txt")[,1]
  XY_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/HumanMethylationSites_X_Y.txt")[,1]
  polymorphic_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/polymorphic_cpgs.txt")[,1]
  exclude <- union(nonspecific_probes,union(XY_probes,polymorphic_probes))
  
  C.refactor <- cov.data
  X.refactor <- X[setdiff(rownames(X),exclude),]
  refactor.cov.mdl <-  TCA::refactor(X=X.refactor, k=length(cell.types),
                                     C=C.refactor, verbose=verbose, 
                                     rand_svd = TRUE,
                                     log_file=NULL)
  refactor.sites <- rownames(refactor.cov.mdl$coeffs)
  intersecting.sites <- intersect(refactor.sites, available.ref.sites)
  message(sprintf("%i sites intersecting with %i ref sites", length(intersecting.sites), length(available.ref.sites)))
 
  # Run TCA models
  # TCA with sites chosen by refactor
  tca.default.mdl <- tca(X=X[refactor.sites,], W=epidish.est.props, C1=C1, C2=cov.data,
                         refit_W=TRUE, refit_W.features=refactor.sites,
                         refit_W.sparsity=length(refactor.sites),
                         verbose=verbose, log_file=NULL, vars.mle=FALSE,
                         constrain_mu=TRUE)
  
  # TCA with reference sites
  tca.ref.features.mdl <- tca(X=X[available.ref.sites,], W=epidish.est.props, C1=C1, C2=cov.data, 
                              refit_W=TRUE, refit_W.features=available.ref.sites, 
                              refit_W.sparsity=length(available.ref.sites), verbose=verbose,
                              log_file=NULL, vars.mle=FALSE, constrain_mu=TRUE)
  
  # TCA with sites chosen by refactor and noisy initialization
  tca.default.noisy.mdl <- tca(X=X[refactor.sites,], W=noisy.epidish.props, C1=C1, C2=cov.data,
                               refit_W=TRUE, refit_W.features=refactor.sites,
                               refit_W.sparsity=length(refactor.sites),
                               verbose=verbose, log_file=NULL, vars.mle=FALSE,
                               constrain_mu=TRUE)
  
  # TCA with reference sites and noisy initialization
  tca.ref.features.noisy.mdl <- tca(X=X[available.ref.sites,], W=noisy.epidish.props, C1=C1, C2=cov.data, 
                                    refit_W=TRUE, refit_W.features=available.ref.sites, verbose=verbose,
                                    refit_W.sparsity=length(available.ref.sites),
                                    log_file=NULL, vars.mle=FALSE, constrain_mu=TRUE)
  
  rmse <- function(x, y){
    indices <- !(is.na(x) | is.na(y))
    rmse.val <- sqrt(sum((y[indices] - x[indices])^2)/length(indices))
    return(rmse.val)
  }
  
  prop.estimates <- list("EpiDISH"=epidish.est.props,
                         "TCA1 Ref"=tca.ref.features.mdl$W,
                         "TCA1"=tca.default.mdl$W,
                         "EpiDISH2"=noisy.epidish.props,
                         "TCA2 Ref"=tca.ref.features.noisy.mdl$W,
                         "TCA2"=tca.default.noisy.mdl$W)
  
  res <- data.frame()
  for (method in names(prop.estimates)){
    for (cell.type in cell.types){
      prop.cor <- cor(prop.estimates[[method]][samples,cell.type], 
                      W.facs[samples,cell.type], 
                      use="complete.obs")
      prop.rmse <- rmse(prop.estimates[[method]][samples,cell.type], 
                        W.facs[samples,cell.type])
      new.row <- data.frame("Method"=method, "Cell Type"=cell.type, 
                            "PCC"=prop.cor, "RMSE"=prop.rmse)
      res <- rbind(res, new.row)
    }
  }
  
  return(list(prop.estimates=prop.estimates, 
              performance.metrics=res))
}

#' Plot performance metrics for BBC and Koestler datasets
#' 
#' @param bbc.results List output by compare_cell_fraction_estimates() for BBC dataset
#' @param koestler.results List output by compare_cell_fraction_estimates() for Koestler dataset
#' @param bbc.sample.size Samples used in BBC dataset
#' @param koestler.sample.size Samples used in Koestler dataset
plot_results <- function(bbc.results, koestler.results, 
                         bbc.sample.size, koestler.sample.size){
  label.angle <- 45
  bbc.pcc <- ggplot2::ggplot(data=bbc.results$performance.metrics, ggplot2::aes(x=Method, y=PCC, fill=Method)) +
    ggplot2::geom_bar(stat="identity", color="black", 
                      width = 0.7, position = ggplot2::position_dodge(width = 2),
                      show.legend=FALSE) +
    ggplot2::facet_wrap(~Cell.Type, nrow=1, ncol=6) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = label.angle, hjust = 1),
                   axis.title.x = ggplot2::element_blank())
  bbc.rmse <- ggplot2::ggplot(data=bbc.results$performance.metrics, ggplot2::aes(x=Method, y=RMSE, fill=Method)) +
    ggplot2::geom_bar(stat="identity", color="black", 
                      width = 0.7, position = ggplot2::position_dodge(width = 2),
                      show.legend=FALSE) +
    ggplot2::facet_wrap(~Cell.Type, nrow=1, ncol=6) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = label.angle, hjust = 1),
                   axis.title.x = ggplot2::element_blank())
  koestler.pcc <- ggplot2::ggplot(data=koestler.results$performance.metrics, ggplot2::aes(x=Method, y=PCC, fill=Method)) +
    ggplot2::geom_bar(stat="identity", color="black", 
                      width = 0.7, position = ggplot2::position_dodge(width = 2),
                      show.legend=FALSE) +
    ggplot2::facet_wrap(~Cell.Type, nrow=1, ncol=6) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = label.angle, hjust = 1),
                   axis.title.x=element_blank())
  koestler.rmse <- ggplot2::ggplot(data=koestler.results$performance.metrics, ggplot2::aes(x=Method, y=RMSE, fill=Method)) +
    ggplot2::geom_bar(stat="identity", color="black", 
                      width = 0.7, position = ggplot2::position_dodge(width = 2),
                      show.legend=FALSE) +
    ggplot2::facet_wrap(~Cell.Type, nrow=1, ncol=6) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = label.angle, hjust = 1),
                   axis.title.x=element_blank())
  
  title1 <- ggdraw() + 
    draw_label(sprintf("Koestler Data\nSample Size = %i", koestler.sample.size),
               x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  title2 <- ggdraw() + 
    draw_label(sprintf("BBC Data\nSample Size = %i", bbc.sample.size),
               x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  top_row <- plot_grid(koestler.pcc, koestler.rmse,
                       nrow = 2, ncol = 1)
  bottom_row <- plot_grid(bbc.pcc, bbc.rmse,
                          nrow = 2, ncol = 1)
  
  refit.w.plot <- plot_grid(title1, top_row, title2, bottom_row, 
                            label_size = 12, ncol = 1,
                            rel_heights = c(0.1, 1, 0.1, 1))
  return(refit.w.plot)
}

#' preprocessing for BBC dataset
#' 
#' @param idat_dir Directory storing IDAT files
#' @param data_dir Directory to store output beta values.
#' @return Path to the beta value csv
epic_qc <- function(idat_dir, data_dir){
  out.file <- sprintf("%s/BBC_beta_ENmix_QC.csv", data_dir)
  rgSet<-readidat(idat.dir,
                  recursive=T, verbose=T)
  #neg: will use 600 chip internal controls probes to estimate background distribution parameters.
  #RELIC: REgression on Logarithm of Internal Control probes
  background = preprocessENmix(rgSet, bgParaEst="neg", dyeCorr="RELIC",nCores=1)
  ##quantile3: will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together
  normalized = norm.quantile(background, method="quantile3")
  
  #Probe design type bias correction using Regression on Correlated Probes (RCP) method
  probecorrected = rcp(normalized)
  
  ##get betas
  bg_beta = data.frame(probecorrected)
  write.csv(bg_beta, out.file,
            quote=FALSE) 
}

#' Replicate cell type fraction analysis of Koestler et al. and BBC datasets
#'
#' @param koestler List generated by get_koestler_data() or provided in CellTypeSpecificMethylationData
#'                 package under the same name
#' @param bbc List provided in CellTypeSpecificMethylationData package under the same name
#' @param results_dir Directory to save plot of results
#' @param plot_type Extension for saving plot graphics, such as pdf or png
#' @param random_seed Integer seed for replication
#' @return A list of results. Slot \code{results} contains a list of proportions and metrics generated
#'         by each method on each dataset. Slot \code{plot} contains the plot object
#' @export
refit_w_comparison <- function(koestler, bbc, results_dir, plot_type, random_seed=1000){
  koestler.results <- compare_cell_fraction_estimates(X = koestler$X, 
                                                      W.facs = koestler$W.facs,
                                                      C1 = NULL, 
                                                      C2 = koestler$C2[,"plate",drop=FALSE],
                                                      random_seed=random_seed)
  
  bbc.results <- compare_cell_fraction_estimates(X = bbc$X,
                                                 W.facs = bbc$W.facs,
                                                 C1 = NULL,
                                                 C2 = bbc$C2[,"Plate",drop=FALSE],
                                                 random_seed=random_seed)
  
  refit.w.plot <- plot_results(bbc.results, koestler.results, 
                               bbc.sample.size=ncol(bbc$X), 
                               koestler.sample.size=ncol(koestler$X))
  outfile <- paste(results_dir, "/refit_w_comparison.", plot_type,
                   sep="")
  ggplot2::ggsave(outfile, plot = refit.w.plot,
                  width = 10, height = 10, units = "in")
  return(list(koestler.results=koestler.results,
              bbc.results=bbc.results,
              refit.w.plot=refit.w.plot))
}
