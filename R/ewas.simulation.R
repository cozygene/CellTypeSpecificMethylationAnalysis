
run_tca.simulation <- function(bulk){
  n <- ncol(bulk$X)
  m <- nrow(bulk$X)
  k <- ncol(bulk$W)
  
  ## (1)  Run association testing by modeling cell-type-specific effects of the phenotype on the methylation levels.
  
  tca.qvals = matrix(1, m, k)
  rownames(tca.qvals) <- rownames(bulk$X)
  if (bulk$model.direction == 1){
    tca.mdl <- tca(bulk$X,bulk$W, C1 = bulk$y, verbose = FALSE)
    tca.estimates <- tca.mdl$gammas_hat
    tca.pvals <- tca.mdl$gammas_hat_pvals
  }else{
    tca.estimates <- matrix(0,m,k)
    tca.pvals <- matrix(0,m,k)
    for (j in 1:m){
      tca.mdl <- tca(bulk$X[j,,drop=F],bulk$W, C1 = bulk$y[,j,drop=F], verbose = FALSE)
      tca.pvals[j,] <- tca.mdl$gammas_hat_pvals
      tca.estimates[j,] <- tca.mdl$gammas_hat
    }
  }
  for (h in 1:k) tca.qvals[,h] <- p.adjust(tca.pvals[,h], method = "fdr")
  
  ## (2) Run association testing by modeling effect of cell-type-specific methylation on the phenotype (i.e. using tcareg)
  tca.mdl <- tca(bulk$X,bulk$W, verbose = FALSE)
  # marginal test (i.e. testing for each cell type separately while *not* accounting for the effects of other cell types)
  tcareg1.pvals = matrix(1, m, k)
  tcareg1.estimates = matrix(0, m, k)
  tcareg1.qvals = matrix(1, m, k)
  rownames(tcareg1.estimates) <- rownames(bulk$X)
  rownames(tcareg1.qvals) <- rownames(bulk$X)
  rownames(tcareg1.qvals) <- rownames(bulk$X)
  if (bulk$model.direction == 1){
    tcareg1.res <- tcareg(bulk$X, tca.mdl, bulk$y, test = "marginal", save_results = FALSE, log_file = NULL, verbose = FALSE)
    for (h in 1:k) tcareg1.estimates[,h] <- tcareg1.res[[h]]$beta
    for (h in 1:k) tcareg1.pvals[,h] <- tcareg1.res[[h]]$pvals
    for (h in 1:k) tcareg1.qvals[,h] <- p.adjust(tcareg1.res[[h]]$pvals, method = "fdr")
  }else{
    for (j in 1:m){
      tcareg1.res <- tcareg(bulk$X[j,,drop=F], tcasub(tca.mdl, features = rownames(bulk$X)[j], verbose = FALSE), bulk$y[,j], C3 = bulk$W[,1:(k-1)], test = "marginal",save_results = FALSE, log_file = NULL, verbose = FALSE)
      for (h in 1:k) tcareg1.estimates[j,h] <- tcareg1.res[[h]]$beta
      for (h in 1:k) tcareg1.pvals[j,h] <- tcareg1.res[[h]]$pvals
    }
    for (h in 1:k) tcareg1.qvals[,h] <- p.adjust(tcareg1.pvals[,h], method = "fdr")
  }
  # conditional marginal test (i.e. testing for each cell type separately while accounting for the effects of other cell types)
  tcareg2.pvals = matrix(1, m, k)
  tcareg2.qvals = matrix(1, m, k)
  tcareg2.estimates = matrix(0, m, k)
  rownames(tcareg2.pvals) <- rownames(bulk$X)
  rownames(tcareg2.qvals) <- rownames(bulk$X)
  rownames(tcareg2.estimates) <- rownames(bulk$X)
  if (bulk$model.direction == 1){
    tcareg2.res <- tcareg(bulk$X, tca.mdl, bulk$y, test = "marginal_conditional",save_results = FALSE, log_file = NULL, verbose = FALSE)
    tcareg2.estimates <- tcareg2.res$beta
    tcareg2.pvals <- tcareg2.res$pvals
    tcareg2.qvals <- tcareg2.res$qvals
  }else{
    for (j in 1:m){
      tcareg2.res <- tcareg(bulk$X[j,,drop=F], tcasub(tca.mdl, features = rownames(bulk$X)[j], verbose = FALSE), bulk$y[,j], C3 = bulk$W[,1:(k-1)], test = "marginal_conditional",save_results = FALSE, log_file = NULL, verbose = FALSE)
      tcareg2.estimates[j,] <- tcareg2.res$beta
      tcareg2.pvals[j,] <- tcareg2.res$pvals
      tcareg2.estimates[j,] <- tcareg2.res$beta
    }
    for (h in 1:k) tcareg2.qvals[,h] <- p.adjust(tcareg2.pvals[,h], method = "fdr")
  }

  return(list("tca.qvals" = tca.qvals, "tca.pvals" = tca.pvals, "tca.estimates" = tca.estimates, "tcareg1.qvals" = tcareg1.qvals, "tcareg1.pvals" = tcareg1.pvals, "tcareg1.estimates" = tcareg1.estimates, "tcareg2.qvals" = tcareg2.qvals, "tcareg2.pvals" = tcareg2.pvals, "tcareg2.estimates" = tcareg2.estimates))
  
}


run_celldmc.simulation <- function(bulk){
  names <- colnames(bulk$W)
  if (bulk$model.direction == 1){
    res <- CellDMC(bulk$X, bulk$y, bulk$W, adjPMethod = "fdr")
    pvals <- matrix(1,nrow(bulk$X),ncol(bulk$W))
    qvals <- matrix(1,nrow(bulk$X),ncol(bulk$W))
    estimates <- matrix(1,nrow(bulk$X),ncol(bulk$W))
    for (i in 1:ncol(bulk$W)){
      pvals[,i] <- res$coe[[names[i]]]$p
      qvals[,i] <- res$coe[[names[i]]]$adjP
      estimates[,i] <- res$coe[[names[i]]]$Estimate
    }
  }else{
    m <- nrow(bulk$X)
    k <- ncol(bulk$W)
    pvals <- matrix(1,m,k)
    qvals <- matrix(1,m,k)
    estimates <- matrix(1,m,k)
    for (j in 1:m){
      res <- CellDMC(bulk$X[c(j,j),,drop=F], bulk$y[,j,drop=F], bulk$W, adjPMethod = "fdr") # run twice on the j-th site to address a bug in celldmc that doesnt allow to run it on a single site
      for (i in 1:ncol(bulk$W)){
        pvals[j,i] <- res$coe[[names[i]]]$p[1]
        estimates[j,i] <- res$coe[[names[i]]]$Estimate[1]
      }
    }
    for (h in 1:k) qvals[,h] <- p.adjust(pvals[,h], method = "fdr")
  }
  return(list("pvals" = pvals,"qvals" = qvals, "estimates" = estimates))
}


