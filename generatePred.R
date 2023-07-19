#'
#' @description generatePred() fits models using Elastic Net(glmnet), Penalized Elastic Net S-estimator (pense), 
#'              Sparse least trimmed squares regression (robustHD), Robust Stepwise Split Regression (robStepSplitReg),
#'              Robust Stepwise Split Regression with an ensemble, SR Lars, SR Lars with an ensemble
#'              and Huber Loss Regression(hqreg) and stores the MSPE, recall, precision and computing time for each.
#' 
#' @param simdata Output from generateData().
#' @param n_models Number of models for ensemble. 

# Required source files
source("RC_PR.R")
source("SparseShootingS/sparseShootingS.R")

generatePred <- function(simdata, n_models, ...) {
  
  # Size of training data set
  N <- length(simdata$training_data$xtrain) 
  # Number of active parameters
  p.active <- simdata$pactive 
  # Sample size for one training set
  n <- simdata$n
  # Total number of predictors 
  p <- simdata$p
  # Test data
  xtestdata <- simdata$testing_data$xtest 
  ytestdata <- simdata$testing_data$ytest
  
  # Initialization of output array
  pred_output <- array(dim = c(14, 4, N))
  colnames(pred_output) <- c("MSPE","RC", "PR", "CPU")
  rownames(pred_output) <- c("ElasticNet", "DDC_ElasticNet", 
                             "sparseShooting",
                             "Pense", "sparseLTS",
                             "HuberEN",
                             "robStepSplitReg", "robStepSplitRegSelect", 
                             "robStepSplitRegEnsemble", "robStepSplitRegEnsembleSelect", 
                             "srlars", "srlarsSelect",
                             "srlarsEnsemble", "srlarsEnsembleSelect")

  
  for(i in 1:N) {
    
    # Elastic net
    en_final <- tryCatch({
      
      en_cpu <- system.time(en_output <- glmnet::cv.glmnet(x = simdata$training_data$xtrain[[i]],
                                     y = simdata$training_data$ytrain[[i]],
                                     alpha = 3/4))
      MSPE_en <- mean((predict(en_output,xtestdata, lambda = "lambda.min")- ytestdata)^2)/simdata$sigma^2
      coef_en <- coef(en_output,lambda = "lambda.min")[-1]
      RC_en <- RC_PR(coef_en, simdata$active_ind)$rc
      PR_en <- RC_PR(coef_en, simdata$active_ind)$pr
      CPU_en <- en_cpu["elapsed"]
      
      c(MSPE_en, RC_en, PR_en, CPU_en)
      
    }, error = function(e){
      return(c(NA, NA, NA, NA))
    })
    pred_output["ElasticNet",, i] <- en_final
    
    # PENSE
    pense_final <- tryCatch({
      
      pense_cpu <- system.time(pense_output <- pense::adapense_cv(x = simdata$training_data$xtrain[[i]],
                                                                  y=simdata$training_data$ytrain[[i]],
                                                                  alpha = 0.75, cv_k = 10, cv_repl = 1, ...))
      MSPE_pense <- mean((predict(pense_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
      coef_pense <- coef(pense_output)[-1]
      RC_pense <- RC_PR(coef_pense, simdata$active_ind)$rc
      PR_pense <- RC_PR(coef_pense, simdata$active_ind)$pr
      CPU_pense <- pense_cpu["elapsed"]
      
      c(MSPE_pense, RC_pense, PR_pense, CPU_pense)
    }, error = function(e){
      return(c(NA, NA, NA, NA))
    })
    pred_output["Pense",, i] <- pense_final
    
    # Huber-EN
    huber_final <- tryCatch({
      
      HUBER_cpu <- system.time(huber_output <- hqreg::cv.hqreg(simdata$training_data$xtrain[[i]],
                                                               simdata$training_data$ytrain[[i]],
                                                               alpha = 0.75, method = "huber"))
      MPSE_huber <- mean((predict(huber_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
      coef_huber <- coef(huber_output)[-1]
      RC_huber <- RC_PR(coef_huber, simdata$active_ind)$rc
      PR_huber <- RC_PR(coef_huber, simdata$active_ind)$pr
      CPU_huber  <- HUBER_cpu["elapsed"]
      
      c(MPSE_huber, RC_huber, PR_huber, CPU_huber)
      
    }, error = function(e){
      return(c(NA, NA, NA, NA))
    })
    pred_output["HuberEN",, i] <- huber_final
    
    # sparse LTS
    sparseLTS_final <- tryCatch({
      
      lambda_max <- robustHD::lambda0(simdata$training_data$xtrain[[i]], simdata$training_data$ytrain[[i]])
      lambda_grid = rev(exp(seq(log(1e-2*lambda_max),log(lambda_max),length = 50)))
      sparseLTS_cpu <- system.time(sparseLTS_output <- robustHD::sparseLTS(x = simdata$training_data$xtrain[[i]], 
                                                                           y = c(simdata$training_data$ytrain[[i]]), 
                                                                           lambda = lambda_grid,
                                                                           mode = "lambda",
                                                                           tol = 1e-2))
      MSPE_sparseLTS <- mean((predict(sparseLTS_output, xtestdata) - ytestdata)^2)/simdata$sigma^2
      coef_sparseLTS <- coef(sparseLTS_output)[-1]
      RC_sparseLTS <- RC_PR(coef_sparseLTS, simdata$active_ind)$rc
      PR_sparseLTS <- RC_PR(coef_sparseLTS, simdata$active_ind)$pr
      CPU_sparseLTS <- sparseLTS_cpu["elapsed"]
      
      c(MSPE_sparseLTS, RC_sparseLTS, PR_sparseLTS, CPU_sparseLTS)
      
    }, error = function(e){
      return(c(NA, NA, NA, NA))
    }) 
    pred_output["sparseLTS",, i] <- sparseLTS_final
    
    # DDC Elastic Net
    DDC_ElasticNet_final <- tryCatch({
      
      DDCxy <- cellWise::DDC(cbind(simdata$training_data$xtrain[[i]],simdata$training_data$ytrain[[i]]), 
                             DDCpars = list(fastDDC = TRUE, silent = TRUE))
      x_imp <- DDCxy$Ximp[, -ncol(DDCxy$Ximp)]
      y_imp <- DDCxy$Ximp[, ncol(DDCxy$Ximp)]
      
      DDC_en_cpu <- system.time(enDDC_output <- glmnet::cv.glmnet(x = x_imp,
                                                                  y = y_imp,
                                                                  alpha = 0.75))
      MSPE_DDC_ElasticNet<- mean((predict(enDDC_output,xtestdata, lambda = "lambda.min")- ytestdata)^2)/simdata$sigma^2
      coef_DDC_ElasticNet <- coef(enDDC_output,lambda = "lambda.min")[-1]
      RC_DDC_ElasticNet <- RC_PR(coef_DDC_ElasticNet, simdata$active_ind)$rc
      PR_DDC_ElasticNet <- RC_PR(coef_DDC_ElasticNet, simdata$active_ind)$pr
      CPU_DDC_ElasticNet <- DDC_en_cpu["elapsed"]
      
      c(MSPE_DDC_ElasticNet, RC_DDC_ElasticNet, PR_DDC_ElasticNet, CPU_DDC_ElasticNet)
      
    }, error = function(e){
      return(c(NA, NA, NA, NA))
    })
    pred_output["DDC_ElasticNet",, i] <- DDC_ElasticNet_final
    
    # SparseShootingS
    sparseshootingS_final <- tryCatch({
      
      sparseS_cpu <- system.time(sparseS_output <- sparseshooting(x =  simdata$training_data$xtrain[[i]], 
                                                                  y =  simdata$training_data$ytrain[[i]],  
                                                                  wvalue = 3, nlambda = 100))
      sparseS_preds <- sparseS_output$coef[1] + xtestdata %*% sparseS_output$coef[-1]
      MSPE_sparseS <- mean((sparseS_preds -ytestdata)^2)/simdata$sigma^2
      coef_sparseS <- sparseS_output$coef[-1]
      RC_sparseS <- RC_PR(coef_sparseS, simdata$active_ind)$rc
      PR_sparseS <- RC_PR(coef_sparseS, simdata$active_ind)$pr
      CPU_sparseS <- sparseS_cpu["elapsed"]
      
      c(MSPE_sparseS, RC_sparseS, PR_sparseS, CPU_sparseS)
      
    }, error = function(e){
      return(c(NA, NA, NA, NA))
    })
    pred_output["sparseShooting",, i] <- sparseshootingS_final

    # robStepSplitReg
    robStepSplitReg_final <- tryCatch({
      
      robStepSplitReg_cpu <- system.time(robStepSplitReg_output <- robStepSplitReg::robStepSplitReg(x = simdata$training_data$xtrain[[i]], 
                                                                 y=simdata$training_data$ytrain[[i]],
                                                                 model_saturation = c("fixed","p-value")[1],
                                                                 alpha = 0.05,
                                                                 model_size = n-1,
                                                                 compute_coef = TRUE))
                                                      
      MPSE_robStepSplitReg <- mean((predict(robStepSplitReg_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
      coef_robStepSplitReg <- coef(robStepSplitReg_output)[-1]
      RC_robStepSplitReg <- RC_PR(coef_robStepSplitReg, simdata$active_ind)$rc
      PR_robStepSplitReg <- RC_PR(coef_robStepSplitReg, simdata$active_ind)$pr
      CPU_robStepSplitReg <- robStepSplitReg_cpu["elapsed"]
     
      c(MPSE_robStepSplitReg, RC_robStepSplitReg, PR_robStepSplitReg,CPU_robStepSplitReg)
      
    }, error = function(e){
      return(c(NA, NA, NA,NA))
        
    })
    pred_output["robStepSplitReg",, i] <- robStepSplitReg_final
    
    # robStepSplitReg:  variables selection RC & PR check
    robStepSplitRegSelect_final <- tryCatch({
      
      MSPE_robStepSplitRegSelect <- NA
      RC_robStepSplitRegSelect <- sum(unlist(robStepSplitReg_output$selections) %in% simdata$active_ind)/length(simdata$active_ind)
      PR_robStepSplitRegSelect <- sum(unlist(robStepSplitReg_output$selections) %in% simdata$active_ind)/length(unlist(robStepSplitReg_output$selections))
      
      c(MSPE_robStepSplitRegSelect, RC_robStepSplitRegSelect, PR_robStepSplitRegSelect, CPU_robStepSplitReg)
      
      }, error = function(e){
        return(c(NA, NA, NA, NA))
        
    })
    pred_output["robStepSplitRegSelect",, i] <- robStepSplitRegSelect_final
    
    
    # robStepSplitReg Ensemble
    robStepSplitReg_ensemble_final <- tryCatch({
      
      robStepSplitReg_ensemble_cpu <- system.time(robStepSplitReg_ensemble_output <- robStepSplitReg::robStepSplitReg(x = simdata$training_data$xtrain[[i]], 
                                                          y=simdata$training_data$ytrain[[i]],
                                                          model_saturation = c("fixed","p-value")[1],
                                                          n_models = n_models,
                                                          alpha = 0.05,
                                                          model_size = min(n-1, floor(p/n_models)),
                                                          compute_coef = TRUE))
      MPSE_robStepSplitReg_ensemble <- mean((predict(robStepSplitReg_ensemble_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
      coef_robStepSplitReg_ensemble <- coef(robStepSplitReg_ensemble_output)[-1]
      RC_robStepSplitReg_ensemble <- RC_PR(coef_robStepSplitReg_ensemble, simdata$active_ind)$rc
      PR_robStepSplitReg_ensemble <- RC_PR(coef_robStepSplitReg_ensemble, simdata$active_ind)$pr
      CPU_robStepSplitReg_ensemble <- robStepSplitReg_ensemble_cpu["elapsed"]
      
      c(MPSE_robStepSplitReg_ensemble, RC_robStepSplitReg_ensemble, PR_robStepSplitReg_ensemble, CPU_robStepSplitReg_ensemble)
      
      }, error = function(e){
        return(c(NA, NA, NA, NA))
      })
    pred_output["robStepSplitRegEnsemble",, i] <- robStepSplitReg_ensemble_final
    
    # robStepSplitReg Ensemble: variables selection RC & PR check
    robStepSplitReg_ensembleSelect_final <- tryCatch({
      
     MPSE_ensembleSelect <- NA
     RC_ensembleSelect_final <- sum(unlist(robStepSplitReg_ensemble_output$selections) %in% simdata$active_ind)/length(simdata$active_ind)
     PR_ensembleSelect_final <- sum(unlist(robStepSplitReg_ensemble_output$selections) %in% simdata$active_ind)/length(unlist(robStepSplitReg_ensemble_output$selections))
      
     c(MPSE_ensembleSelect, RC_ensembleSelect_final, PR_ensembleSelect_final,CPU_robStepSplitReg_ensemble)
     
     }, error = function(e){
      return(c(NA, NA, NA, NA))
    })
    pred_output["robStepSplitRegEnsembleSelect",, i] <- robStepSplitReg_ensembleSelect_final

    # Split Robust LARS
    srlars_final <- tryCatch({
      
      srlars_cpu <- system.time(srlars_output <- srlars::srlars(simdata$training_data$xtrain[[i]],
                          simdata$training_data$ytrain[[i]], 
                          n_models = 1,
                          model_saturation = c("fixed", "p-value")[1],
                          alpha = 0.05, model_size = n-1,
                          robust = TRUE,
                          compute_coef = TRUE,
                          en_alpha = 1/4))
    
      MSPE_srlars<- mean((predict(srlars_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
      coef_srlars <- coef(srlars_output)[-1]
      RC_srlars <- RC_PR(coef_srlars, simdata$active_ind)$rc
      PR_srlars <- RC_PR(coef_srlars, simdata$active_ind)$pr
      CPU_srlars  <- srlars_cpu["elapsed"]
      
      c(MSPE_srlars, RC_srlars, PR_srlars,CPU_srlars)
    
    }, error = function(e){
      return(c(NA, NA, NA, NA))
    })
    pred_output["srlars",, i] <- srlars_final
    
    # Split Robust LARS: variables selection RC & PR check
    srlarsSelect_final <- tryCatch({
        
      MSPE_srlarsSelect <- NA
      RC_srlarsSelect <- sum(unlist(srlars_output$selections) %in% simdata$active_ind)/length(simdata$active_ind)
      PR_srlarsSelect <- sum(unlist(srlars_output$selections) %in% simdata$active_ind)/length(unlist(srlars_output$selections))
        
      c(MSPE_srlarsSelect, RC_srlarsSelect, PR_srlarsSelect,CPU_srlars)
    
    }, error = function(e){
        return(c(NA, NA, NA, NA))
    })
    pred_output["srlarsSelect",, i] <- srlarsSelect_final
  
    
    # Split Robust LARS Ensemble
    srlarsEnsemble_final <- tryCatch({
        
      srlarsEnsemble_cpu <- system.time(srlars_outputEnsemble <- srlars::srlars(simdata$training_data$xtrain[[i]],
                                      simdata$training_data$ytrain[[i]], 
                                      n_models = n_models,
                                      model_saturation = c("fixed", "p-value")[1],
                                      alpha = 0.05, model_size = min(n-1, floor(p/n_models)),
                                      robust = TRUE,
                                      compute_coef = TRUE,
                                      en_alpha = 1/4))
      
      MSPE_srlarsEnsemble<- mean((predict(srlars_outputEnsemble,newx = xtestdata)- ytestdata)^2)/simdata$sigma^2
      coef_srlarsEnsemble <- coef(srlars_outputEnsemble)[-1]
      RC_srlarsEnsemble <- RC_PR(coef_srlarsEnsemble, simdata$active_ind)$rc
      PR_srlarsEnsemble <- RC_PR(coef_srlarsEnsemble, simdata$active_ind)$pr
      CPU_srlarsEnsemble  <- srlarsEnsemble_cpu["elapsed"]
      
      c(MSPE_srlarsEnsemble, RC_srlarsEnsemble, PR_srlarsEnsemble, CPU_srlarsEnsemble)
      
    }, error = function(e){
      return(c(NA, NA, NA, NA))
    })
    pred_output["srlarsEnsemble",, i] <- srlarsEnsemble_final
    
    # Split Robust LARS Ensemble: variable selection PR & RC check 
    srlarsEnsembleSelect_final <- tryCatch({
        
      MSPE_srlarsEnsembleSelect <- NA
      RC_srlarsEnsembleSelect <- sum(unlist(srlars_outputEnsemble$selections) %in% simdata$active_ind)/length(simdata$active_ind)
      PR_srlarsEnsembleSelect <- sum(unlist(srlars_outputEnsemble$selections) %in% simdata$active_ind)/length(unlist(srlars_outputEnsemble$selections))
      
      c(MSPE_srlarsEnsembleSelect, RC_srlarsEnsembleSelect, PR_srlarsEnsembleSelect, CPU_srlarsEnsemble)
      
    }, error = function(e){
      return(c(NA, NA, NA, NA))
    })
    pred_output["srlarsEnsembleSelect",, i] <- srlarsEnsembleSelect_final
  }
  
  # Return full prediction output
  return (pred_output)
}

