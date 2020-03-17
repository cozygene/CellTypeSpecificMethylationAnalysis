library(TCA)
library(ggplot2)

#' Compares sigma and mu esimtation of TCA with or without vars.mle option enabled
compare_vars_mle <- function(sample.sizes, n.replications, n.sites, n.cell.types,
                             n.c1.covs, n.c2.covs, tau, random_seed=1000){
  res <- data.frame()
  set.seed(random_seed)
  for (n.samples in sample.sizes){
    for (replicate in 1:n.replications){
      sim.data <- TCA::test_data(n = n.samples, m = n.sites, k = n.cell.types,
                                 p1 = n.c1.covs, p2 = n.c2.covs, tau = tau, 
                                 log_file = NULL, verbose = FALSE)
      
      fast.mdl <- TCA::tca(X=sim.data$X, W=sim.data$W, 
                           C1=sim.data$C1, C2=sim.data$C2,
                           vars.mle = FALSE,
                           log_file = NULL, verbose = FALSE)
      
      orig.mdl <- TCA::tca(X=sim.data$X, W=sim.data$W, 
                           C1=sim.data$C1, C2=sim.data$C2,
                           vars.mle = TRUE,
                           log_file = NULL, verbose = FALSE)
      
      fast.mus.cor <- sapply(1:n.cell.types, function(h){
        cor(fast.mdl$mus_hat[,h], sim.data$mus[,h])
      })
      names(fast.mus.cor) <- paste0("fast.mus.cor.", 1:n.cell.types)
      fast.sigmas.cor <- sapply(1:n.cell.types, function(h){
        cor(fast.mdl$sigmas_hat[,h], sim.data$sigmas[,h])   
      })
      names(fast.sigmas.cor) <- paste0("fast.sigmas.cor.", 1:n.cell.types)
      orig.mus.cor <- sapply(1:n.cell.types, function(h){
        cor(orig.mdl$mus_hat[,h], sim.data$mus[,h])
      })
      names(orig.mus.cor) <- paste0("orig.mus.cor.", 1:n.cell.types)
      orig.sigmas.cor <- sapply(1:n.cell.types, function(h){
        cor(orig.mdl$sigmas_hat[,h], sim.data$sigmas[,h])   
      })
      names(orig.sigmas.cor) <- paste0("orig.sigmas.cor.", 1:n.cell.types)
      experiment.cols <- c(n.samples, replicate)
      names(experiment.cols) <- c("n.samples", "replicate")
      res.row <- t(data.frame(c(experiment.cols, fast.mus.cor, orig.mus.cor, fast.sigmas.cor, orig.sigmas.cor)))
      res <- rbind(res, res.row)
      message(sprintf("Replicate %i with %i samples... done.", replicate, n.samples))
    }
  }
  rownames(res) <- 1:nrow(res)
  return(res)
}

#' Summarizes result dataframe for ggplot2 barplots
summarize_res <- function(res){
  sample.sizes <- sort(unique(res$n.samples))
  fast.mus.col <- "fast.mus.cor.1"
  fast.sigmas.col <- "fast.sigmas.cor.1"
  orig.mus.col <- "orig.mus.cor.1"
  orig.sigmas.col <- "orig.sigmas.cor.1"
  res_summary <- data.frame()
  for (sample.size in sample.sizes){
    sub_res <- res[res$n.samples == sample.size,]
    fast.mu.mean <- mean(sub_res[,fast.mus.col])
    fast.mu.sd <- sd(sub_res[,fast.mus.col])
    fast.sigma.mean <- mean(sub_res[,fast.sigmas.col])
    fast.sigma.sd <- sd(sub_res[,fast.sigmas.col])
    fast.mu.row <- data.frame(n.samples=sample.size, method="GMM",
                              param="Mean",
                              mean=fast.mu.mean, sd=fast.mu.sd)
    fast.sigma.row <- data.frame(n.samples=sample.size, method="GMM",
                                 param="Variance",
                                 mean=fast.sigma.mean, sd=fast.sigma.sd)
    orig.mu.mean <- mean(sub_res[,orig.mus.col])
    orig.mu.sd <- sd(sub_res[,orig.mus.col])
    orig.sigma.mean <- mean(sub_res[,orig.sigmas.col])
    orig.sigma.sd <- sd(sub_res[,orig.sigmas.col])
    orig.mu.row <- data.frame(n.samples=sample.size, method="MLE",
                              param="Mean",
                              mean=orig.mu.mean, sd=orig.mu.sd)
    orig.sigma.row <- data.frame(n.samples=sample.size, method="MLE",
                                 param="Variance",
                                 mean=orig.sigma.mean, sd=orig.sigma.sd)
    res_summary <- rbind(res_summary, fast.mu.row, fast.sigma.row, orig.mu.row, orig.sigma.row)
  }
  res_summary$n.samples <- factor(as.character(res_summary$n.samples), levels = as.character(sort(unique(res_summary$n.samples))))
  res_summary$param <- factor(res_summary$param, levels= c("Variance", "Mean"))
  res_summary$method <- factor(res_summary$method, levels=c("GMM", "MLE"))
  return(res_summary)
}

#' Generates plot with ggplot2
plot_cors <- function(res.summary){
  facet.labels <- paste("Estimated",c("variance", "mean"), "parameters")
  names(facet.labels) <- c("Variance", "Mean")
  runtimes.plot <- ggplot2::ggplot(data=res.summary, ggplot2::aes(x=n.samples, y=mean, fill=method)) +
                   ggplot2::geom_bar(stat="identity", color="black", position=ggplot2::position_dodge()) +
                   ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-sd, ymax=mean+sd), width=0.1, position=ggplot2::position_dodge(0.9)) +
                   ggplot2::facet_wrap(~param, nrow=1, ncol=2,
                                       labeller=ggplot2::labeller(param=facet.labels)) +
                   ggplot2::theme_bw() + ggplot2::scale_fill_grey(start = .7, end = .9) + 
                   ggplot2::xlab("Sample size") +
                   ggplot2::ylab("PCC")
  return(runtimes.plot)
}

#' Compare parameter estimation of TCA in MLE versus GMM Mode
#'
#' @param results_dir Directory to save plot of results
#' @param plot_type Extension for saving plot graphics, such as pdf or png
#' @return A list of the results and plot
#' @export
var_mle_gmm_comparison <- function(results_dir, plot_type){
  res <- compare_vars_mle(sample.sizes=c(100,500,1000),
                          n.replications=5,
                          n.sites=100,
                          n.cell.types=5,
                          n.c1.covs=2,
                          n.c2.covs=2,
                          tau=0.01,
                          random_seed=1000)

  res.summary <- summarize_res(res)
  res.plot <- plot_cors(res.summary)
  outfile <- paste(results_dir, "/var_mle_gmm_comparison.", plot_type,
                   sep="")
  ggplot2::ggsave(outfile, plot = res.plot, width = 7, height=4.6, units = "in")
  return(list(mle.gmm.comparison.results=res.summary,
              mle.gmm.comparison.plot=res.plot))
}
