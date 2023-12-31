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
  
  if(contamination.scenario %in% c("cellwise_marginal", "cellwise_correlation", "casewise")){
    
    # Creating list for output
    output <- lapply(1:length(contamination.prop), function(t1) return(lapply(1:length(p.active), function(t2) return(list()))))
    
    # Computation of simfunc() for different values of p.active and contamination proportion 
    for(contamination.id in 1:length(contamination.prop)){
      for(active.id in 1:length(p.active)) {
        
        # Print p.active and sparsity
        cat("\n", "tau: ", contamination.prop[contamination.id])
        cat("\n", "p.active: ", p.active[active.id], "\n")
        
        output[[contamination.id]][[active.id]] <- simfunc(N = N, n = n, m = m, p = p, 
                                                           rho = rho, rho.inactive = rho.inactive, 
                                                           p.active = p.active[active.id], 
                                                           group.size = group.size, snr = snr, 
                                                           contamination.prop = contamination.prop[contamination.id], 
                                                           contamination.scenario = contamination.scenario,
                                                           n_models = n_models, 
                                                           ...)
      }
    }
    
  } else if(contamination.scenario %in% c("mixture_marginal", "mixture_correlation")){
    
    # Creating list for output
    output <- lapply(1:length(p.active), function(t1) return(list()))
    
    # Computation of simfunc() for different values of p.active and contamination proportion 
    for(active.id in 1:length(p.active)) {
      
      # Print p.active
      cat("\n", "p.active: ", p.active[active.id], "\n")
      
      output[[active.id]] <- simfunc(N = N, n = n, m = m, p = p, 
                                     rho = rho, rho.inactive = rho.inactive, 
                                     p.active = p.active[active.id], 
                                     group.size = group.size, snr = snr, 
                                     contamination.prop = contamination.prop, 
                                     contamination.scenario = contamination.scenario,
                                     n_models = n_models, 
                                     ...)
      }
  }
  
  # Returning the list of outputs
  return(output)
}

