#simfunc() calls generateData() and generatePred() and produces a matrix containing 

source("generateData.R")
source("generatePred.R")


simfunc <- function(N, #number of training sets
                    n, #sample size for training set
                    m, #sample size for test set
                    parameters, #number of parameters
                    rho, #correlation
                    rho.inactive,
                    p.active, #number of active parameters
                    group.size,
                    snr, 
                    contamination.prop,
                    n_models,
                    ...){ #contamination proportion
  
  sim_Data <- generateData(N, n, m, parameters, rho, rho.inactive, p.active, group.size, snr, contamination.prop)
  output <- generatePred(sim_Data, n_models, ...)
  
  return(output)
  }
