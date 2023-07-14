#' 
#' @description simfunc() calls generateData() and generatePred().
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

# Required source files 
source("generateData.R")
source("generatePred.R")

simfunc <- function(N, 
                    n, 
                    m, 
                    parameters, 
                    rho, 
                    rho.inactive,
                    p.active, 
                    group.size,
                    snr, 
                    contamination.prop,
                    contamination.scenario,
                    n_models,
                    ...){ 
  
  sim_Data <- generateData(N, n, m, parameters, rho, rho.inactive, p.active, group.size, snr, 
                           contamination.prop, contamination.scenario)
  output <- generatePred(sim_Data, n_models, ...)
  
  return(output)
}
