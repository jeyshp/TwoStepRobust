# TwoStepRobust
Simulations scripts for fast and scalable cellwise robust ensembles for high-dimensional data.

#### Required library
```
install.packages("mvnfast")
install.packages("parallel")
install.packages("foreach")
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
contamination_scenario <- c("casewise", "cellwise_marginal", "cellwise_correlation")
snr <- c(0.5, 1, 2)
rho <- c(0.5, 0.8)
n_models <- 2
p.active <- c(50, 100, 200)
```

#### Loading Source Files

`source("generateOutput.R")`

#### Generating Results

```
for(scenario_val in contamination_scenario) {
  
  if(scenario_val == "casewise") 
    contamination.prop = seq(0, 0.4, by = 0.1) else 
      contamination.prop = c(0.1, 0.2)
    
  for(snr_val in snr){
    for(rho_val in rho) {
      filename = paste0("results/results_n=",n,"_p=",p,"_scenario=", scenario_val, "_snr=",snr_val,"_rho= ",rho_val,".Rdata")
      result <- generateOutput(N=N, n=n, m= 2e5, p=p, rho=rho_val, rho.inactive = 0.2,
                               p.active = p.active, group.size= 2, snr = snr_val , 
                               contamination.prop = contamination.prop, contamination.scenario = scenario_val,
                               seed = 0, n_models = n_models)
      save.image(filename)
    
    }
  }
}
```
