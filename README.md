# TwoStepRobust
Simulations scripts for fast and scalable cellwise robust ensembles for high-dimensional data.

#### Required library
```
install.packages("mvnfast")
install.packages("parallel")
install.packages("doSNOW")
install.packages("pense")
install.packages("robustHD")
install.packages("robStepSplitReg")
install.packages("srlars")
install.packages("hqreg")
install.packages("glmnet")
install.packages("pcaPP") 
install.packages("robustbase") 
```

#### Simulation Parameters

```
N <- 50
n <- 50
p <- 500
contamination_scenario <- c("casewise", 
                            "cellwise_marginal", "cellwise_correlation",
                            "mixture_marginal", "mixture_correlation")
snr <- c(0.5, 1, 2)
rho <- 0.8
n_models <- 5
p.active <- c(50, 100, 200)
```

#### Loading Source Files

`source("generateOutput.R")`

#### Generating Results

```
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
```
