# TwoStepRobust
Simulations to test and compare the robStepSplitReg estimator's performance on high-dimensional contaminated data. 

#### Required library
```
library("mvnfast")
library("parallel")
library("doSNOW")
library("pense")
library("robustHD")
library("robStepSplitReg")
library("hqreg")
library("glmnet")
```


#### Simulation Parameters

```
N <- 2
n <- 25
p <- 50
snr = c(0.5,1,2)
rho = c(0.5,0.8)
n_models = 2
```

`source("generateOutput.R")`

#### Generating & saving results for one simulation 

```
for(snr_val in snr){
  for(rho_val in rho) {
    filename = paste0("results/results_n=",n,"_p=",p,"_snr=",snr_val,"_rho= ",rho_val,".Rdata")
    result <- generateOutput(N=N, n=n, m= 2e5, p=p, rho=rho_val, rho.inactive = 0.2,
                             p.active = c(5,10), group.size= 2, snr = snr_val , contamination.prop = c(0.2,0.3),
                             seed = 1464, n_models = n_models)
    save.image(filename)
  
  }
}
```
