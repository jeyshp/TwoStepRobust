#' @description generateOutput() calls simfunc() in parallel for varying values of p.active and contamination prop.
#' 
#' @param N Number of training sets.
#' @param n Sample size for training set.
#' @param m Sample size for test set.
#' @param p Total number of parameters.
#' @param rho Correlation within a block of active parameters.
#' @param rho.inactive Correlation between blocks of active parameters.
#' @param p.active Number of active parameters.
#' @param group.size Size of one block of active parameters.
#' @param snr Signal to noise ratio.
#' @param contamination.prop Contamination proportion.
#' @param n_models Number of models for ensemble.
#' @param contamination.scenario Casewise, cellwise marginal or cellwise correlation contamination.

# Required source file
source("simfunc.R")

generateOutput <- function (N, 
                            n , 
                            m, 
                            p, 
                            rho, 
                            rho.inactive = 0.2,
                            p.active, 
                            group.size,
                            snr, 
                            contamination.prop, 
                            contamination.scenario,
                            seed = 0,
                            n_models,
                            ...){
  # Setting seed
  set.seed(seed)
  
  # Initializing number of cores used 
  if (length(p.active)*length(contamination.prop) != 1) 
    n_clusters <-  min(length(p.active)*length(contamination.prop), detectCores()-1) else
      n_clusters <- 2
  mycluster <- makeCluster(n_clusters)
  registerDoSNOW(mycluster)
  
  if(contamination.scenario %in% c("casewise", "cellwise_marginal", "cellwise_correlation")){
    
    # Parallel computation of simfunc() for different values of p.active and contamination proportion 
    results <- foreach(mycontprop = contamination.prop) %:%
      foreach(mypactive = p.active) %dopar% {
        
        output <- simfunc(N = N, n = n, m = m, p = p, rho = rho, rho.inactive = rho.inactive, p.active = mypactive, 
                          group.size = group.size, snr = snr, 
                          contamination.prop = mycontprop, contamination.scenario = contamination.scenario,
                          n_models = n_models, 
                          ...)
      }
  } else if(contamination.scenario %in% c("mixture_marginal", "mixture_correlation")){
    
    # Parallel computation of simfunc() for different values of p.active and contamination proportion 
    results <- foreach(mypactive = p.active) %dopar% {
        
        output <- simfunc(N = N, n = n, m = m, p = p, rho = rho, rho.inactive = rho.inactive, p.active = mypactive, 
                          group.size = group.size, snr = snr, 
                          contamination.prop = contamination.prop, contamination.scenario = contamination.scenario,
                          n_models = n_models, 
                          ...)
      }
  }
  
  stopCluster(mycluster)
  return(results)
}

