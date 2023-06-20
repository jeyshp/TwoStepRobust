
library("parallel")
install.packages("doSNOW")
library("doSNOW")

#generateOutput() calls simfunc in parallel for varying values of p.active, snr, contamination prop and rho

source("simfunc.R")

generateOutput <- function (N, #number of training sets
                            n , #sample size for training set
                            m, #sample size for test set
                            p, #number of parameters
                            rho, #correlation
                            rho.inactive = 0.2,
                            p.active, #number of active parameters
                            group.size,
                            snr, 
                            contamination.prop, #contamination proportion
                            seed = 0,
                            n_models,
                            ...){
  set.seed(seed)
  if (length(p.active)*length(contamination.prop) != 1) {
    n_clusters =  min(length(p.active)*length(contamination.prop),detectCores()-1)
  } else
    n_clusters = 2
  
  mycluster = makeCluster(n_clusters)
  
  #present results as table (fixed(n.p), fixed(rho,snr)) - varying prop of active p and cont.prop
  #saving files 
  results = foreach(mycontprop = contamination.prop, .packages = c("pense","robustHD","mvnfast","robStepSplitReg","hqreg","glmnet")) %:%
    foreach(mypactive = p.active, .packages = c("pense","robustHD","mvnfast","robStepSplitReg","hqreg","glmnet")) %dopar% {
      output = simfunc(N = N, n=n, m=m, p=p, rho = rho, rho.inactive = rho.inactive, p.active = mypactive, 
                       group.size = group.size, snr = snr, contamination.prop = mycontprop, n_models = n_models, ...)
    }
  
  stopCluster(mycluster)
  return(results)
}

