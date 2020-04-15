
get_reference_data <- function(data_path){
  gse.file <- paste(data_path,"GSE35069_Matrix_signal_intensities.txt",sep="/")
  gse.gz.file <- paste(data_path,"GSE35069_Matrix_signal_intensities.txt.gz",sep="/")
  if (!file.exists(gse.file)){
    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE35069&format=file&file=GSE35069%5FMatrix%5Fsignal%5Fintensities%2Etxt%2Egz", paste(data_path,"GSE35069_Matrix_signal_intensities.txt.gz",sep="/"))
    GEOquery::gunzip(gse.gz.file)
  }
  cell_types.names <- c("Granulocytes", "CD4+", "CD8+", "monocytes","CD19+","NK")
  cell_types.sample_pos <- cbind(seq(13,18),seq(19,24),seq(25,30),seq(31,36),seq(37,42),seq(43,48))
  Z.raw <- read.table(gse.file, sep="\t", header=TRUE, row.names = 1)
  k <- length(cell_types.names)
  Z.beta <- vector(mode="list", length=k)
  for (h in 1:k){
    Z.beta[[h]] <- matrix(0,nrow(Z.raw),6)
    rownames(Z.beta[[h]]) <- rownames(Z.raw)
    counter <- 1
    for (j in cell_types.sample_pos[,h]){
      Z.beta[[h]][,counter] <- Z.raw[,(j-1)*3+2] / (Z.raw[,(j-1)*3+1] + Z.raw[,(j-1)*3+2] +100)
      counter <- counter + 1
    }
  }
  file.remove(gse.file)
  return (Z.beta)  
}


get_hannum <- function(data_path){
  gse.file <- paste(data_path,"GSE40279_average_beta.txt",sep="/")
  if (!file.exists(gse.file)){
    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40279&format=file&file=GSE40279%5Faverage%5Fbeta%2Etxt%2Egz",paste(data_path,"GSE40279_average_beta.txt.gz",sep="/"))
    GEOquery::gunzip(paste(data_path, "GSE40279_average_beta.txt.gz",sep="/"))
  }
  X <- data.table::fread(gse.file,header = TRUE)
  X.cpgs <- as.matrix(X[,1])[,1]
  X <- as.matrix(X[,2:ncol(X)])
  rownames(X) <- X.cpgs
  file.remove(gse.file)
  return (X)
}


init_summary_lists <- function(effect_sizes, num_sims){
  celldmc.summary <- list()
  for (i in 1:4){
    scenario <- paste("scenario",i,sep="")
    celldmc.summary[[scenario]] <- list()
    celldmc.summary[[scenario]][["SE"]] <- matrix(0,length(effect_sizes),num_sims)
    celldmc.summary[[scenario]][["SP"]] <- matrix(0,length(effect_sizes),num_sims)
    celldmc.summary[[scenario]][["PPV"]] <- matrix(0,length(effect_sizes),num_sims)
  }
  tca.methods <- c("tca","tcareg1","tcareg2")
  tca.summary <- list()
  for (method in tca.methods) tca.summary[[method]] <- data.table::copy(celldmc.summary)
  return(list("tca.summary" = tca.summary, "celldmc.summary" = celldmc.summary))
}


get_scenario_effect_pairs <- function(effect_sizes, num_scenarios){
  S <- matrix(0,num_scenarios*length(effect_sizes),2)
  counter <- 1
  for (i in 1:num_scenarios){
    for (j in 1:length(effect_sizes)){
      S[counter,1] <- i
      S[counter,2] <- effect_sizes[j]
      counter <- counter + 1
    }
  }
  return(S)
}

prep_parametric_simulation_data <- function(data_path){
  
  file_name <- paste(data_path, "/parametric_simulation_data.RData", sep="")
  if (!file.exists(file_name)){
    # Download and prepare the Reinius et al. data from which we will estimate the parameters for data simnulation.
    Z.beta <- get_reference_data(data_path)
    
    # Load the reference data that will be used for estimating cell-type proportions
    # blood.ref.333cpgs.txt is the list of cpgs and their means as described in in Teschendorff et al. 2017, BMC Bioinformatics.
    #ref.tmp <- read.table("assets/blood.ref.333cpgs.txt", header = T, row.names = 1)
    #ref.tmp <- EpiDISH::centDHSbloodDMC.m
    #ref <- as.matrix(ref.tmp[,c("Gran","CD4T","CD8T","Mono","B","NK")])
    ref <- blood.ref.333cpgs
    ref.cpgnames <- rownames(ref)
    
    # Download the Hannum et al. data
    X.hannum <- get_hannum(data_path)
    # Estimate cell-type proportions for the samples in the Hannum data; these estimateswill be used as a pool of cell-type proportions for simulating bulk data.
    W.hannum <- epidish(X.hannum, ref)$estF
    rm(X.hannum)
    
    # Estimate the parameters of the refrence data; these will be used for simulating data;
    Z.params <- fit_beta(Z.beta)
    Z.s1 <- Z.params$Z.s1
    Z.s2 <- Z.params$Z.s2
    Z.mean <- Z.params$Z.mean
    Z.var <- Z.params$Z.var
    
    # Dont allow variances and means to be exactly zero (to avoid precision issues)
    Z.var[Z.var == 0] <- 1e-7
    Z.mean[Z.mean == 0] <- 1e-7
    
    save(Z.beta, ref, ref.cpgnames, W.hannum, Z.s1, Z.s2, Z.mean, Z.var, file = file_name)
  }
  
  return(file_name)
}


parametric.run_simulation <- function(n, m, m.true, effect_sizes, num_sims, num_cores, model.direction, require_effect_direction, parametric_simulation_data_file, results_path, random_seed=1000){
  set.seed(random_seed)
  
  load(parametric_simulation_data_file)
  
  k <- ncol(W.hannum)
  non_ref.cpgnames <- setdiff(rownames(Z.s1),ref.cpgnames)
  
  lst <- init_summary_lists(effect_sizes, num_sims)
  celldmc.summary <- lst$celldmc.summary
  tca.summary <- lst$tca.summary
  
  tryCatch({
    
    cl <- NULL
    if (num_cores-1){
      cl <- makeCluster(num_cores)
      invisible(clusterEvalQ(cl, c(library("nloptr"),library("testit"),library("pracma"),library("TCA"),library("EpiDISH") )))  
    } 
    
    for (num_sim in 1:num_sims){
      sim_success <- FALSE
      tries <- 0
      max.tries <- 3
      while (!sim_success){
        print(num_sim)
        # simulate cell-type-specific levels; add the 333 reference sites at the end
        Z.sim <- parametric.simulate_Z(Z.beta, Z.s1, Z.s2, ref.cpgnames, non_ref.cpgnames, n, m) 
        Z.sim.cpgnames <- rownames(Z.sim[[1]])
        # scenario-effect size pairs to run
        S <- get_scenario_effect_pairs(effect_sizes, 4)
        
        if (num_cores-1) clusterExport(cl, varlist = c("k","m","S","Z.sim","Z.s1","Z.s2","Z.mean","Z.var",
                                                       "Z.sim.cpgnames","W.hannum","parametric.simulate_bulk",
                                                       "run_celldmc.simulation","run_tca.simulation","evaluate_celldmc_results",
                                                       "evaluate_tca_results", "parametric.select_true_sites", "calculate_metrics",
                                                       "m.true","model.direction","ref", "require_effect_direction"), envir=environment())
        # ADDED tryCatch here for CellDMC bug when feature has too little variance
        tryCatch({ 
          res <- pbapply::pblapply(1:nrow(S),function(j){
            scenario <- S[j,1];
            effect_size <- S[j,2];
            bulk <- parametric.simulate_bulk(Z.sim, Z.s1[Z.sim.cpgnames,], Z.s2[Z.sim.cpgnames,], Z.mean[Z.sim.cpgnames,], Z.var[Z.sim.cpgnames,], ref, W.hannum, m.true, effect_size, as.character(scenario), model.direction = model.direction);
            celldmc.res <- run_celldmc.simulation(bulk);
            # if an error occurs at celldmc.res, restart this loop (at most max.tries times)
            celldmc.metrics <- evaluate_celldmc_results(celldmc.res, bulk$true.cpgs, m, k, require_effect_direction);
            tca.res <- run_tca.simulation(bulk);
            tca.metrics <- evaluate_tca_results(tca.res, bulk$true.cpgs, m, k, require_effect_direction);
            res.metrics <- list();
            res.metrics[["celldmc.res"]] <- celldmc.metrics;
            res.metrics[["tca.res"]] <- tca.metrics;
            return(res.metrics) },cl = cl )
          sim_success <- TRUE},
          error = function(c){message("There was an error!")
                              message(c)})
        if (!sim_success){
          tries <- tries + 1
          if (tries >= max.tries){
            stop("Maximum attempts for simulation reached, simulation parameters may be causing issues with CellDMC.")
          }
        }
      }
      lst <- update_metrics_summaries(S, res, num_sim, tca.summary, celldmc.summary)
      tca.summary <- lst$tca.summary
      celldmc.summary <- lst$celldmc.summary
    }
    
    #if (parametric_simulation_data_file){
    if (require_effect_direction){
      filename <- paste(results_path,"/parametric_simulation_results_m_",m,"_n_",n,"_model_direction_",model.direction,".require_effect_direction.RData",sep="")
    }else{
      filename <- paste(results_path,"/parametric_simulation_results_m_",m,"_n_",n,"_model_direction_",model.direction,".RData",sep="")
    }
    save(celldmc.summary, tca.summary, effect_sizes, n, m, m.true, model.direction, file = filename)
    
  }, finally = {
    if (num_cores-1) stopCluster(cl)
  })
  
}


# simulate cell-type-specific levels
parametric.simulate_Z <- function(Z.beta, Z.s1, Z.s2, ref.cpgnames, non_ref.cpgnames, n, m){
  k <- ncol(Z.s1)
  cpgs <- c(non_ref.cpgnames[sample(length(non_ref.cpgnames))[1:m]],ref.cpgnames)
  Z <- vector(mode="list", length=k)
  for (h in 1:k){
    Z[[h]] <- matrix(0,length(cpgs),n)
    rownames(Z[[h]]) <- cpgs
    for (i in 1:length(cpgs)){
      Z[[h]][i,] <- rbeta(n,Z.s1[cpgs[i],h],Z.s2[cpgs[i],h])
    }
  }
  return(Z)
}

fit_beta <- function(Z.beta){
  k <- length(Z.beta)
  cpgnames <- rownames(Z.beta[[1]])
  m <- nrow(Z.beta[[1]])
  Z.mean <- matrix(0,m,k)
  rownames(Z.mean) <- cpgnames
  Z.var <- matrix(0,m,k)
  rownames(Z.var) <- cpgnames
  Z.s1 <- matrix(0,m,k)
  rownames(Z.s1) <- cpgnames
  Z.s2 <- matrix(0,m,k)
  rownames(Z.s2) <- cpgnames
  for (h in 1:k){
    for (j in 1:m){
      if (sd(Z.beta[[h]][j,]) == 0) next
      res <- EnvStats::ebeta(Z.beta[[h]][j,], method = "mme")  
      Z.s1[j,h] <- res$parameters[1]
      Z.s2[j,h] <- res$parameters[2]
      s <- Z.s1[j,h] + Z.s2[j,h]
      Z.mean[j,h] <- Z.s1[j,h]/s
      Z.var[j,h] <- Z.mean[j,h]*(1-Z.mean[j,h])/(s+1)
    }
  }
  return(list("Z.s1" = Z.s1, "Z.s2" = Z.s2, "Z.mean" = Z.mean, "Z.var" = Z.var))
}


parametric.select_true_sites <- function(Z.mean, Z.var, m.true, m.true.cell_types, effect_size, mode, ref.cpgnames, exclude.ref){
  m.true.low = 0.2
  m.true.high = 0.8
  opts <- list("algorithm"="NLOPT_GN_DIRECT","xtol_rel"=1.0e-8)
  if (effect_size == 0) m.true <- 0
  true.cpgs <- list()
  true.cpgs[["cpg"]] <- numeric(m.true)
  true.cpgs[["effect_size"]] <- numeric(m.true) + effect_size
  true.cpgs[["cell_types"]] <- matrix(0,m.true,m.true.cell_types)
  true.cpgs[["changes"]] <- matrix(0,m.true,m.true.cell_types)
  true.cpgs[["deltas"]] <- matrix(0,m.true,m.true.cell_types)
  true.cpgs[["beta.s1"]] <- matrix(0,m.true,m.true.cell_types)
  true.cpgs[["beta.s2"]] <- matrix(0,m.true,m.true.cell_types)
  if (effect_size == 0) return (true.cpgs)
  # select candidate sites
  if (mode == "unidirectional"){
    true.cpgs.pool <- which( rowSums(Z.mean < m.true.low) >= m.true.cell_types |  rowSums(Z.mean > m.true.high) >= m.true.cell_types )
  }else{
    true.cpgs.pool <- which(rowSums(Z.mean < m.true.low) > 0 & rowSums(Z.mean > m.true.high) > 0 & rowSums(Z.mean < m.true.low | Z.mean > m.true.high) >= m.true.cell_types)
  }
  
  # select true sites
  perm <- sample(length(true.cpgs.pool))
  for (i in 1:m.true){
    cpg <- perm[i]
    true.cpgs$position[i] <- true.cpgs.pool[cpg]
    pool.cell_types <- which(Z.mean[true.cpgs.pool[cpg],] < m.true.low | Z.mean[true.cpgs.pool[cpg],] > m.true.high)
    if (mode == "unidirectional"){
      while(sum(true.cpgs$cell_types[i,]) == 0){
        cell_types <- pool.cell_types[sample(length(pool.cell_types))[1:m.true.cell_types]]
        means <- Z.mean[true.cpgs.pool[cpg],cell_types]
        if (sum(means < m.true.low) ==  m.true.cell_types | sum(means > m.true.high) == m.true.cell_types){
          true.cpgs$cell_types[i,] <- cell_types
        }
      }
    }else{
      while(sum(true.cpgs$cell_types[i,]) == 0){
        cell_types <- pool.cell_types[sample(length(pool.cell_types))[1:m.true.cell_types]]
        means <- Z.mean[true.cpgs.pool[cpg],cell_types]
        if (sum(means < m.true.low) > 0 & sum(means > m.true.high) > 0){
          true.cpgs$cell_types[i,] <- cell_types
        }
      }
    }
    # parameters of the controls group
    mu2 <- Z.mean[true.cpgs.pool[cpg],true.cpgs$cell_types[i,]]
    sigma2 <- (Z.var[true.cpgs.pool[cpg],true.cpgs$cell_types[i,]])**0.5
    for (h in 1:length(mu2)){
      # solve a non-linear equation with a single parameter in order to get delta; constrain mu1 to be hypo or hyper methylated based on mu2
      func <- function(x) (abs(mu2[h]-x)/sqrt((sigma2[h]**2 + (x*(1-x)/(mu2[h]*(1-mu2[h])))*sigma2[h]**2)/2) - effect_size )**2
      if (mu2[h] < m.true.low){
        true.cpgs$changes[i,h] <- 1 # 1 for hyper
        mu1 <- nloptr::nloptr(mu2[h], func, lb = mu2[h], ub = 1, opts = opts)$solution
        testit::assert(mu1 > mu2[h])
      }else{
        true.cpgs$changes[i,h] <- -1 # -1 for hypo
        mu1 <- nloptr::nloptr(mu2[h], func, lb = 0, ub = mu2[h], opts = opts)$solution
        testit::assert(mu1 < mu2[h])
      }
      sigma1 <- sqrt(mu1*(1-mu1)/(mu2[h]*(1-mu2[h])))*sigma2[h]
      true.cpgs$beta.s1[i,h] <- ((1-mu1)/(sigma1**2) - 1/mu1)*(mu1**2)
      true.cpgs$beta.s2[i,h] <- true.cpgs$beta.s1[i,h]*(1/mu1-1)
      true.cpgs$deltas[i,h] <- abs(mu1-mu2[h])
    }
  }
  return(true.cpgs)
}


parametric.simulate_bulk <- function(Z, Z.s1, Z.s2, Z.mean, Z.var, ref, W.pool, m.true, effect_size, true.mode, model.direction){
  n <- ncol(Z[[1]])
  n.controls <- n/2
  n.cases <- n- n.controls
  m <- nrow(Z[[1]])
  
  # Sample cell-type proportions
  W <- W.pool[sample(nrow(W.pool), n, replace = T),]
  k <- ncol(W)
  
  ref.cpgnames <- rownames(ref)
  if (true.mode == "1") true.cpgs <- parametric.select_true_sites(Z.mean, Z.var, m.true, 1, effect_size, mode = "unidirectional", ref.cpgnames)
  if (true.mode == "2") true.cpgs <- parametric.select_true_sites(Z.mean, Z.var, m.true, 2, effect_size, mode = "unidirectional", ref.cpgnames)
  if (true.mode == "3") true.cpgs <- parametric.select_true_sites(Z.mean, Z.var, m.true, 2, effect_size, mode = "bidirectional", ref.cpgnames)
  if (true.mode == "4") true.cpgs <- parametric.select_true_sites(Z.mean, Z.var, m.true, 3, effect_size, mode = "bidirectional", ref.cpgnames)
  
  if (model.direction){
    # cell-type-specific methylation is affected by the phenotype
    y <- matrix(0,n,1)
    rownames(y) <- rownames(W)
    colnames(y) <- "1"
    y[(n.controls+1):n] <- 1
    for (j in 1:m.true){
      cell_types <- true.cpgs$cell_types[j,]
      for (h in 1:length(cell_types)){
        Z[[true.cpgs$cell_types[j,h]]][true.cpgs$position[j],(n.controls+1):n] <- rbeta(n.cases,true.cpgs$beta.s1[j],true.cpgs$beta.s2[j])  
      }
    }
  }else{
    # For each true site generate a phenotype that is affected by the site; note that in principle we could have generated a single phenotype that is affected by multiple sites, however, one phenotype per associated site allows a more clear interpretation of the effect sizes.
    y <- matrix(0,n,m)
    colnames(y) <- 1:m
    rownames(y) <- rownames(W)
    for (j in 1:m.true){
      # set the variance of y to be the sum of variances of the cell types for which we assume effect sizes.
      total_var <- 0
      for (cell_type in true.cpgs$cell_types[j,]) total_var <- total_var + Z.var[true.cpgs$position[j],cell_type]
      y[,true.cpgs$position[j]] <- rnorm(n, mean = 0, sd = sqrt(total_var)) + W%*%rnorm(ncol(W), mean = 0, sd = 1)
      for (h in 1:length(true.cpgs$cell_types[j,])){
        y[,true.cpgs$position[j]] <- y[,true.cpgs$position[j]] + effect_size * Z[[true.cpgs$cell_types[j,h]]][true.cpgs$position[j],] * true.cpgs$changes[j,h]
      }
    }
    # for each non-true site randomly pick one of the phenotypes for the association test.
    for (j in setdiff(1:m,true.cpgs$position)){
      y[,j] <- y[,true.cpgs$position[sample(m.true)[1]]]
    }
  }
  
  X <- matrix(0,m,n)
  for (h in 1:k) X <- X + Z[[h]]*pracma::repmat(W[,h],m,1)
  colnames(X) <- rownames(W)
  
  # Estimate W using a reference-based approach
  W.estimated <- epidish(X, ref)$estF
  
  return(list("Z" = Z, "W.real" = W, "W" = W.estimated, "X" = X, "y" = y, "true.cpgs" = true.cpgs, "model.direction" = model.direction))
}

plot_power_simulation <- function(outfile, methods, methods.names, effect_sizes){

  titles <- c(rep(c("Uni-1C"),3),rep(c("Uni-2C"),3),rep(c("Bi-2C"),3),rep(c("Bi-3C"),3))
  metrics <- c("SE","SP","PPV")
  scenarios <- c("scenario1","scenario2","scenario3","scenario4")
  
  plots <- list(mode = vector, length = length(scenarios)*length(metrics))
  counter <- 1
  ylims <- list()
  for (metric in metrics){
    ylims[[counter]] <- c(1,0)
    for (scenario in scenarios){
      for (i in 1:length(methods)){
        min_val <- min(methods[[i]][[scenario]][[metric]])
        if (min_val < ylims[[counter]][1]) ylims[[counter]][1] <- min_val
        max_val <- max(methods[[i]][[scenario]][[metric]])
        if (max_val > ylims[[counter]][2]) ylims[[counter]][2] <- max_val
      }
    }
    counter <- counter + 1
  }
  counter <- 1
  for (scenario in scenarios){
    metric.counter <- 1
    for (metric in metrics){
      num_sims <- ncol(methods[[1]][[scenario]][[metric]])
      data <- data.frame(matrix(ncol = 3, nrow = length(methods)*length(effect_sizes)*num_sims))
      colnames(data) <- c("method","effect_size","value")
      for (i in 1:length(methods)){
        data[((i-1)*length(effect_sizes)*num_sims+1):(i*length(effect_sizes)*num_sims),] <- cbind(rep(c(methods.names[i]),length(effect_sizes)*num_sims), rep(effect_sizes,num_sims), Reshape(methods[[i]][[scenario]][[metric]],length(effect_sizes)*num_sims,1))
      }
      data$effect_size <- as.factor(data$effect_size)
      data$value <- as.numeric(data$value)
      data2 <- data.frame(matrix(ncol = 3, nrow = length(effect_sizes)*length(methods)))
      for (i in 1:length(methods)){
        data2[((i-1)*length(effect_sizes)+1):(i*length(effect_sizes)),] <- cbind(rep(c(methods.names[i]),length(effect_sizes)), effect_sizes,rowMeans(methods[[i]][[scenario]][[metric]]))
      }
      colnames(data2) <- c("method","effect_size","median")
      data2$median <- as.numeric(data2$median)
      data2$effect_size <- as.numeric(data2$effect_size)
      p <- ggplot(data, aes(x=effect_size, y=value, fill=method)) + geom_violin(alpha = 0.25) + theme_ipsum(axis_title_just="c") + xlab("Effect size") + ylab(metric) + ggtitle(titles[counter]) + theme(plot.margin=unit(c(0,0,0.5,0.5), "cm")) + theme(legend.position="bottom") + theme(plot.title = element_text(hjust = 0.5),legend.margin=margin(t = 0, unit='cm'),axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), legend.title=element_blank()) + ylim(ylims[[metric.counter]])
      plots[[counter]] <- p + geom_line(data=data2, aes(x=effect_size, y=median, group=method, color=method),size=1, alpha = 0.8)
      counter <- counter + 1
      metric.counter <- metric.counter + 1
    }
  }
  p.final <- ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],plots[[9]],plots[[10]],plots[[11]],plots[[12]],
            ncol = 3, nrow = 4,
            labels = c("a","b","c","d","e","f","g","h","i","g","k","l"))
  ggsave(outfile,plot = p.final, width = 12, height = 12)
}

#' Replicate parametric simulations for power analyses
#'
#' @param data_dir Path to directory to store data
#' @param results_dir Path to directory to store results
#' @param plot_dir Path to directory to store plots
#' @param plot_type Type of image to save plot
#' @param experiment_index An integer in 1, 2, 3, or 4. 1 specifies phenotype
#'                         affecting cell-type-specific methylation levels.
#'                         2 specifies the same as one but true positives
#'                         also consider predicted direction.
#'                         3 specifies the phenotype being affected by
#'                         cell-type-specific methylation levels. 4 is the
#'                         same as 1 with smaller sample size
#' @param num_cores Number of cores to use for analysis
#' @param m Number CpGs to simulate
#' @param m.true Number of true associations
#' @param effect_sizes Vector of effect sizes to simulate
#' @param n Sample size
#' @param num_sims Number of simulations to run
#' @return None (TODO: Return a list of results)
#' @export
parametric_simulations <- function(data_dir, results_dir, plot_dir,
                                   plot_type = "tiff",
                                   experiment_index=c(1,2,3,4), num_cores=1,
                                   m=1000, m.true=100, effect_sizes=1:8,
                                   n=500, num_sims=50){
  random_seed <- 1000
  if (m.true > m){
    stop("m.true should be less than or equal to m")
  }
  # for consistent variable names across functions
  plots_path <- plot_dir
  results_path <- results_dir
  image_format <- plot_type
  # download and prepare data
  parametric_simulation_data_file <- prep_parametric_simulation_data(data_dir)
  if (experiment_index == 1){
    # run simulations that model the phenotype as affecting methylation at the cell-type level (by setting model.direction = 1)
    parametric.run_simulation(n, m, m.true, effect_sizes, num_sims, num_cores, model.direction = 1, require_effect_direction = FALSE, parametric_simulation_data_file, results_path, random_seed)
    # plots
    load(paste(results_path,"/parametric_simulation_results_m_",m,"_n_",n,"_model_direction_1.RData",sep=""))
    plot_power_simulation(paste(plots_path,"/Figure1.", image_format, sep=""), list(celldmc.summary, tca.summary$tca), c("CellDMC","TCA (X|Y)"), effect_sizes)
    plot_power_simulation(paste(plots_path,"/FigureS7.", image_format, sep=""), list(celldmc.summary, tca.summary$tcareg1), c("CellDMC","TCA (Y|X, marginal)"), effect_sizes)
    
  }
  if (experiment_index == 2){
    # run simulations that model the phenotype as affecting methylation at the cell-type level (by setting model.direction = 1); this time consider a predicted true associations as a true positive only if it captured the true direction of effect.
    parametric.run_simulation(n, m, m.true, effect_sizes, num_sims, num_cores, model.direction = 1, require_effect_direction = TRUE, parametric_simulation_data_file, results_path, random_seed)
    
    # plots
    load(paste(results_path,"/parametric_simulation_results_m_",m,"_n_",n,"_model_direction_1.require_effect_direction.RData",sep=""))
    plot_power_simulation(paste(plots_path,"/FigureS8.", image_format, sep=""), list(celldmc.summary, tca.summary$tcareg1), c("CellDMC","TCA (Y|X, marginal)"), effect_sizes)
  }
  if (experiment_index == 3){
    # run simulations that model the phenotype as affected by cell-type-specific methylation (by setting model.direction = 0)
    parametric.run_simulation(n, m, m.true, effect_sizes, num_sims, num_cores, model.direction = 0, require_effect_direction = FALSE, parametric_simulation_data_file, results_path, random_seed)
    
    # plots
    load(paste(results_path,"/parametric_simulation_results_m_",m,"_n_",n,"_model_direction_0.RData",sep=""))
    plot_power_simulation(paste(plots_path,"/Figure2.",image_format, sep=""), list(celldmc.summary,tca.summary$tcareg2), c("CellDMC","TCA (Y|X)"), effect_sizes)
    plot_power_simulation(paste(plots_path,"/FigureS1.", image_format, sep=""), list(celldmc.summary, tca.summary$tca), c("CellDMC","TCA (X|Y)"), effect_sizes)
    plot_power_simulation(paste(plots_path,"/FigureS9.",image_format, sep=""), list(celldmc.summary,tca.summary$tcareg1), c("CellDMC","TCA (Y|X, marginal)"), effect_sizes)
  }
  if (experiment_index == 4){
    # rerun experiment 1 with smaller sample size
    parametric.run_simulation(30, m, m.true, effect_sizes, num_sims, num_cores, model.direction = 1, require_effect_direction = FALSE, parametric_simulation_data_file, results_path, random_seed)
    load(paste(results_path,"/parametric_simulation_results_m_",m,"_n_30_model_direction_1.RData",sep=""))
    plot_power_simulation(paste(plots_path,"/FigureS3.", image_format, sep=""), list(celldmc.summary, tca.summary$tca), c("CellDMC","TCA (X|Y)"), effect_sizes)
    
    parametric.run_simulation(60, m, m.true, effect_sizes, num_sims, num_cores, model.direction = 1, require_effect_direction = FALSE, parametric_simulation_data_file, results_path, random_seed)
    load(paste(results_path,"/parametric_simulation_results_m_",m,"_n_60_model_direction_1.RData",sep=""))
    plot_power_simulation(paste(plots_path,"/FigureS4.", image_format, sep=""), list(celldmc.summary, tca.summary$tca), c("CellDMC","TCA (X|Y)"), effect_sizes)
  }
}
