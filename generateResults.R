#'
#' @description Simulation example.
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

# Clear all memory
rm(list = ls())

# Required libraries
.libPaths(c("/home/jeysh/R/x86_64-pc-linux-gnu-library/4.3" , "/opt/R/4.3.0/lib/R/library"))

library(parallel)

# Required source file
source("generateOutput.R")

# Setting values for size of training set, sample size and total number of parameters
N <- 50
n <- 50
p <- 500
contamination_scenario <- c("cellwise_marginal", "cellwise_correlation",
                            "mixture_marginal", "mixture_correlation",
                            "casewise")
snr <- c(0.5, 1, 2)
rho <- 0.8
n_models <- 5
p.active <- c(50, 100, 200)

# For loops for varying values of snr and rho a 
for(scenario_val in contamination_scenario) {
  
  if(scenario_val == "casewise") 
    contamination.prop <- seq(0, 0.2, by = 0.1) else if (scenario_val %in% c("cellwise_marginal", 
                                                                             "cellwise_correlation")) 
      contamination.prop <- seq(0.05, 0.1, by = 0.05) else if(scenario_val %in% c("mixture_marginal", 
                                                                                  "mixture_correlation"))
        contamination.prop <- c(0.1, 0.05)
    
  for(snr_val in snr){
    for(rho_val in rho) {
      
      # Print simulation information
      cat("\n", "Contamination Scenario: ", scenario_val)
      cat("\n", "SNR: ", snr_val)
      cat("\n", "Within-Block Correlation: ", rho_val, "\n")
      
      # File name with specifications of simulation
      filename <- paste0("results/results_n=",n,"_p=",p,"_scenario=", scenario_val, 
                         "_snr=", snr_val,"_rho= ", rho_val, ".Rdata")
      
      # Generating results of simulation
      results <- generateOutput(N = N, n = n, m = 2e5, p = p, 
                                rho = rho_val, rho.inactive = 0.2,
                                p.active = p.active, group.size = 25, 
                                snr = snr_val, 
                                contamination.prop = contamination.prop, 
                                contamination.scenario = scenario_val,
                                seed = 0, n_models = n_models)
      
      # Saving simulation results
      save.image(filename)
    }
  }
}

#' @description get_specific() returns the  output for MSPE, RC and PR for all models for a specific p.active ad contamination proportion 
#' 
#' @param result Output from generateOutput()
#' @param p.active Number of active predictors to check for 
#' @param p.active_vec List of active predictors fed to generateOutput()
#' @param contamination.prop Contamination proportion to check for 
#' @param contamination.prop_vec List of contamination proportion to check for 
#' 
get_specific <- function(result,
                         p.active,contamination.prop,
                         p.active_vec,
                         contamination.prop_vec){
  
  return(result[[which(p.active == p.active_vec)]][[which(contamination.prop == contamination.prop_vec)]])
}


