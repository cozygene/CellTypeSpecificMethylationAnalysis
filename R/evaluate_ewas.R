
calculate_metrics <- function(true.cpgs, pvals, estimated_effects,  require_effect_direction,th = 0.05){
  k <- ncol(pvals)
  true.cpgs.rownames <- rownames(true.cpgs)
  cpgs.rownames <- rownames(pvals)
  tps <- 0
  fps <- 0
  # go over truly associated sites
  for (l in 1:length(true.cpgs$position)){
    for (h in 1:k){
      positive <- pvals[true.cpgs$position[l],h] < th
      if (h %in% true.cpgs$cell_types[l,]){
        direction <- true.cpgs$changes[l,which(true.cpgs$cell_types[l,] == h)] == estimated_effects[true.cpgs$position[l],h]/abs(estimated_effects[true.cpgs$position[l],h])
        # whether to require the estimated effect size to match in direction with the true effect size
        if (require_effect_direction){
          if (positive & direction) tps <- tps + 1; if (positive & !direction) fps <- fps + 1  
        }else{
          if (positive) tps <- tps + 1  
        }
      }else{
        if (positive) fps <- fps + 1
      }
    }
  }
  # count fps among the truly non-associated sites
  fps <- fps + sum(pvals[setdiff(1:nrow(pvals),true.cpgs$position),] < th)
  # summarize metrics
  positives <- nrow(true.cpgs$cell_types)*ncol(true.cpgs$cell_types)
  SE <- tps/positives
  SP <- 1-fps/(nrow(pvals)*ncol(pvals) - positives)
  total_positions <- tps+fps
  PPV <- total_positions
  if (total_positions) PPV <- tps/total_positions
  return(list("SE" = SE, "SP" = SP, "PPV" = PPV))
}


evaluate_tca_results <- function(tca.res, true.cpgs, m, k, require_effect_direction){
  res <- list()
  tca.methods <- c("tca","tcareg1","tcareg2")
  for (method in tca.methods){
    res[[method]] <- list()
    l <- calculate_metrics(true.cpgs, tca.res[[paste(method,".qvals",sep="")]], tca.res[[paste(method,".estimates",sep="")]], require_effect_direction)
    res[[method]] <- list()
    res[[method]][["SE"]] <- l$SE
    res[[method]][["SP"]] <- l$SP
    res[[method]][["PPV"]] <- l$PPV
  }
  return(res)
}

evaluate_celldmc_results <- function(celldmc.res, true.cpgs, m, k, require_effect_direction){
  res <- list()
  l <- calculate_metrics(true.cpgs, celldmc.res$qvals, celldmc.res$estimates, require_effect_direction)
  res[["SE"]] <- l$SE
  res[["SP"]] <- l$SP
  res[["PPV"]] <- l$PPV
  return(res)
}


update_metrics_summaries <- function(S, res, num_sim, tca.summary, celldmc.summary){
  
  prev_scenario <- S[1,1]
  effect_size.index <- 0
  for (j in 1:length(res)){
    scenario <- paste("scenario",as.character(S[j,1]),sep="")
    if (S[j,1] != prev_scenario){
      effect_size.index <- 1
      prev_scenario = S[j,1]
    }else{
      effect_size.index <- effect_size.index + 1
    }
    celldmc.summary[[scenario]]$SP[effect_size.index,num_sim] <- res[[j]]$celldmc.res$SP
    celldmc.summary[[scenario]]$SE[effect_size.index,num_sim] <- res[[j]]$celldmc.res$SE
    celldmc.summary[[scenario]]$PPV[effect_size.index,num_sim] <- res[[j]]$celldmc.res$PPV
    for (method in ls(tca.summary)){
      tca.summary[[method]][[scenario]]$SP[effect_size.index,num_sim] <- res[[j]]$tca.res[[method]]$SP
      tca.summary[[method]][[scenario]]$SE[effect_size.index,num_sim] <- res[[j]]$tca.res[[method]]$SE
      tca.summary[[method]][[scenario]]$PPV[effect_size.index,num_sim] <- res[[j]]$tca.res[[method]]$PPV
    }
  }
  
  return(list("tca.summary" = tca.summary, "celldmc.summary" = celldmc.summary))
  
}
