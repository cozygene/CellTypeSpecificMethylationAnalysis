library(TCA)
library(EpiDISH)
library(ggplot2)
library(grid)

#' Run original TCA method
#'
#' methods vector should have these exact function names. Input should always
#' be X, W, Y as produced by the simulation framework
tca.mle <- function(X,W,Y){
         tca.mdl <- TCA::tca(X=X, W=W,
                             vars.mle=TRUE,
                             num_cores=1,
                             log_file=NULL)
         return(0)
       }

#' Run fast TCA method
#'
#' methods vector should have these exact function names. Input should always
#' be X, W, Y as produced by the simulation framework
tca.gmm <- function(X,W,Y){
             tca.mdl <- TCA::tca(X=X, W=W,
                                 vars.mle=FALSE,
                                 num_cores=1,
                                 log_file=NULL)
             return(0)
           }

#' Run CellDMC method
#'
#' methods vector should have these exact function names. Input should always
#' be X, W, Y as produced by the simulation framework
celldmc <- function(X,W,Y){
             EpiDISH::CellDMC(beta.m=X, pheno.v=Y, frac.m=W,
                              mc.cores=1)
             return(0)
           }

#' Compare method runtimes across various input sizes
#'
#' @param sample.sizes Numeric vector for sample sizes to simulate.
#' @param number.sites Numeric vector for number of features/sites to simulate.
#' @param methods Character vector of method wrapper names to benchmark.
#'                Examples of wrappers (tca.mle, tca.gmm, celldmc) are above. They
#'                must have the same input arguments. Function names must
#'                match to be called correctly.
#' @param method.names Pretty names of those listed in methods argument for plotting
#' @param number.replications Integer indicating how many replications of
#'                            each experiment to run. Same input is run for
#'                            each replicate. This is meant to account for
#'                            variability in runtime due to the system we are
#'                            running on.
#' @param number.cell.types Integer indicating number of cell types to simulate.
#' @param seed Integer to set seed for simulation to ensure replicability
#' @return A list of dataframes. First dataframe with columns indicating the method,
#'         replicate number, sample size, number of sites, number of sources, and
#'         total runtime for each experiment. Second dataframe summarizes replicates.
compare_runtimes <- function(sample.sizes, number.sites, methods, method.names,
                             number.replications, number.cell.types, seed){
  runtimes <- data.frame()
  for (n in sort(sample.sizes, decreasing = TRUE)){
    for (m in sort(number.sites, decreasing = TRUE)){
      set.seed(seed)
      simtime <- system.time(sim.data <- TCA::test_data(n, m, number.cell.types, 
                                                        p1=1, p2=0, tau=0.1,
                                                        verbose=F))["elapsed"]
      X <- sim.data[["X"]]
      W <- sim.data[["W"]]
      colnames(W) <- paste0("CELLTYPE", 1:number.cell.types)
      Y <- sim.data[["C1"]]
      message(sprintf("Starting %i samples and %i sites", n, m))
      for (replicate in 1:number.replications){
        for (method.name in methods){
          runtime <- system.time(match.fun(method.name)(X,W,Y))["elapsed"]
          row <- data.frame(method=method.name, replicate=replicate, sample.size=n,
                            number.sites=m, number.sources=number.cell.types,
                            runtime=as.numeric(runtime))
          runtimes <- rbind(runtimes, row)
        }
        message(sprintf("Replicate %i of %i done", 
                        replicate, number.replications))
      }
    }
  }

  # Summarize runtimes over replicates
  runtimes.summary <- data.frame()
  for (method.index in 1:length(methods)){
    method <- methods[method.index]
    method.name <- method.names[method.index]
    for (sample.size in sample.sizes){
      for (num.sites in number.sites){
        runtime.experiment <- runtimes[runtimes$method == method & 
                                         runtimes$sample.size == sample.size & 
                                         runtimes$number.sites == num.sites,]
        runtime.mean <- mean(runtime.experiment$runtime)
        runtime.sd <- sd(runtime.experiment$runtime)
        runtime.experiment <- data.frame(Method=method.name, sample.size=sample.size, 
                                         number.sites=num.sites, runtime.mean=runtime.mean, 
                                         runtime.sd=runtime.sd)
        runtimes.summary <- rbind(runtimes.summary, runtime.experiment)
      }
    }
  }
  return(list(runtimes=runtimes, runtimes.summary=runtimes.summary))
}

#' Plots runtime summary dataframe
plot_runtimes <- function(runtimes.summary){
  sample.sizes <- sort(unique(runtimes.summary[["sample.size"]]))
  number.sites <- sort(unique(runtimes.summary[["number.sites"]]))
  add.error.bar <- TRUE
  if (anyNA(runtimes.summary[["runtime.sd"]])){
    add.error.bar <- FALSE
    y.limits <- c(min(runtimes.summary[["runtime.mean"]]), 
                  max(runtimes.summary[["runtime.mean"]]))
  }
  else {
    y.limits <- c(min(runtimes.summary[["runtime.mean"]] - runtimes.summary[["runtime.sd"]]), 
                  max(runtimes.summary[["runtime.mean"]] + runtimes.summary[["runtime.sd"]]))
  }

  y.breaks <- nchar(as.character(round(y.limits)))
  y.breaks <- 10^(y.breaks[1]:y.breaks[2])
  facet.labels <- paste(sample.sizes, "samples")
  names(facet.labels) <- sample.sizes
  runtimes.plot <- ggplot2::ggplot(runtimes.summary, 
                                   ggplot2::aes(x=number.sites,
                                                y=runtime.mean,
                                                group=Method)) +
    ggplot2::facet_wrap(~sample.size, nrow=1, ncol=length(sample.sizes),
                        labeller=ggplot2::labeller(sample.size=facet.labels)) + 
    ggplot2::geom_point(ggplot2::aes(color=Method)) + 
    ggplot2::geom_line(ggplot2::aes(color=Method)) + 
    ggplot2::scale_x_continuous(trans='log10', 
                                breaks=number.sites, 
                                labels=paste0(number.sites/1000, 'k'), 
                                name="Number of CpGs") +
    ggplot2::scale_y_continuous(trans='log10', 
                                limits=y.limits, 
                                breaks=y.breaks,
                                name="Runtime (seconds)") + 
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="bottom",
                   axis.ticks = ggplot2::element_blank(), 
                   aspect.ratio=1, plot.margin=grid::unit(c(1,1,1,1), "mm"))
  if (add.error.bar){
    runtimes.plot <- runtimes.plot + 
      ggplot2::geom_errorbar(ggplot2::aes(ymin=runtime.mean-runtime.sd,
                             ymax=runtime.mean + runtime.sd,
                             color=Method), width=.05)
  }
  return(runtimes.plot)
}

#' Compare runtimes of CellDMC and TCA (GMM and MLE modes)
#'
#' @param results_dir Directory to save plot of results
#' @param plot_type Extension for saving plot graphics, such as pdf or png
#' @return A list of the runtime results and plot
#' @export
compare_runtimes <- function(results_dir, plot_type){
  runtime.results <- compare_runtimes(sample.sizes = c(50, 100, 200),
                                      number.sites = c(1000, 10000, 50000, 100000),
                                      methods = c("celldmc", "tca.gmm", "tca.mle"),
                                      method.names = c("CellDMC", "TCA (GMM)", "TCA (MLE)"),
                                      number.replications = 1,
                                      number.cell.types = 6,
                                      seed = 42)

  runtime.plot <- plot_runtimes(runtime.results$runtimes.summary)
  outfile <- paste(results_dir, "/compare_runtimes.", plot_type,
                   sep="")
  ggplot2::ggsave(outfile, plot = runtime.plot,  width = 10, height=4.6, units = "in")
  return(list(runtime.results=runtime.results,
              runtime.plot=runtime.plot))
}
