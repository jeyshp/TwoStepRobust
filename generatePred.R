
#' @description generatePred() fits models using Elastic Net(glmnet), Penalized Elastic Net S-estimator (pense), 
#' Sparse least trimmed squares regression (robustHD), Robust Stepwise Split Regression (robStepSplitReg),
#' Robust Stepwise Split Regression with an ensemble and Huber Loss Regression(hqreg) and stores the MSPE, recall and precision for each.
#' 
#' @param simdata Output from generateData() 
#' @param n_models Number of models for ensemble 
#' 


source("RC_PR.R")

generatePred <- function(simdata,n_models, ...) {
  
  #Size of training data set
  N = length(simdata$training_data$xtrain) 
  #Number of active parameters
  p.active = simdata$pactive 
  #Sample size for one training set
  n = simdata$n
  #Total number of predictors 
  p = simdata$p
  #Test data
  xtestdata = simdata$testing_data$xtest 
  ytestdata = simdata$testing_data$ytest
  
  #Initialization of output array
  pred_output <- array(dim = c(8,3,N))
  colnames(pred_output) = c("MSPE","RC", "PR")
  rownames(pred_output) = c("EN", "Pense", "sparseLTS","robStepSplitRegSelect", "robStepSplitReg",
                            "EnsembleSelect","Ensemble", "HuberEN")
  
  
  
  for(i in 1:N) {
    
    #GLMNET
    en_final <- tryCatch({
      en_output <- glmnet::cv.glmnet(x = simdata$training_data$xtrain[[i]],
                                     y = simdata$training_data$ytrain[[i]],
                                     alpha = 0.75)
      MSPE_en <- mean((predict(en_output,xtestdata, lambda = "lambda.min")- ytestdata)^2)/simdata$sigma^2
      coef_en <- coef(en_output,lambda = "lambda.min")[-1]
      RC_en <- RC_PR(coef_en, simdata$active_ind)$rc
      PR_en <- RC_PR(coef_en, simdata$active_ind)$pr
      
      c(MSPE_en, RC_en, PR_en)
      
    }, error = function(e){
      return(c(NA, NA, NA))
    })
    
    pred_output["EN",,i] <- en_final
    
    #PENSE
    pense_final <- tryCatch(
      {
        pense_output <- pense::adapense_cv(x = simdata$training_data$xtrain[[i]],
                                           y=simdata$training_data$ytrain[[i]],
                                           alpha = 0.75, cv_k = 10, cv_repl = 1, ...)
        MSPE_pense <- mean((predict(pense_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
        coef_pense <- coef(pense_output)[-1]
        RC_pense <- RC_PR(coef_pense, simdata$active_ind)$rc
        PR_pense <- RC_PR(coef_pense, simdata$active_ind)$pr
        
        c(MSPE_pense, RC_pense, PR_pense)
      }, error = function(e){
        return(c(NA, NA, NA))
        
     })
    pred_output["Pense",,i] <- pense_final
    
    
  #sparseLTS
    sparseLTS_final <- tryCatch(
      {
        lambda_max <- robustHD::lambda0(simdata$training_data$xtrain[[i]], simdata$training_data$ytrain[[i]])
        lambda_grid = rev(exp(seq(log(1e-2*lambda_max),log(lambda_max),length = 50)))
        sparseLTS_output <- robustHD::sparseLTS(x = simdata$training_data$xtrain[[i]], 
                                                y = c(simdata$training_data$ytrain[[i]]), 
                                                lambda = lambda_grid,
                                                mode = "lambda",
                                                tol = 1e-2)
        MSPE_sparseLTS <- mean((predict(sparseLTS_output, xtestdata) - ytestdata)^2)/simdata$sigma^2
        coef_sparseLTS <- coef(sparseLTS_output)[-1]
        RC_sparseLTS <- RC_PR(coef_sparseLTS, simdata$active_ind)$rc
        PR_sparseLTS <- RC_PR(coef_sparseLTS, simdata$active_ind)$pr
        
        c(MSPE_sparseLTS, RC_sparseLTS, PR_sparseLTS)
      }, error = function(e){
        return(c(NA, NA, NA))
        
      }) 
    pred_output["sparseLTS",,i] <- sparseLTS_final
    
    #robStepSplitReg
    robStepSplitReg_final <- tryCatch(
      {
        robStepSplitReg_output <- robStepSplitReg::robStepSplitReg(x = simdata$training_data$xtrain[[i]], 
                                                                   y=simdata$training_data$ytrain[[i]],
                                                                   model_saturation = c("fixed","p-value")[1],
                                                                   alpha = 0.05,
                                                                   model_size = n-1,
                                                                   compute_coef = TRUE,
                                                                   enpy_opts = pense::enpy_options(retain_max = 50),
                                                                   eps = 1e-1)
        MPSE_robStepSplitReg <- mean((predict(robStepSplitReg_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
        coef_robStepSplitReg <- coef(robStepSplitReg_output)[-1]
        RC_robStepSplitReg <- RC_PR(coef_robStepSplitReg, simdata$active_ind)$rc
        PR_robStepSplitReg <- RC_PR(coef_robStepSplitReg, simdata$active_ind)$pr
        
       
        c(MPSE_robStepSplitReg, RC_robStepSplitReg, PR_robStepSplitReg)
      }, error = function(e){
        return(c(NA, NA, NA))
        
      }
    )
    pred_output["robStepSplitReg",,i] <- robStepSplitReg_final
    
    #robStepSplitReg First variables selection through step wise 
    
    robStepSplitRegSelect_final <- tryCatch(
      {
      
      MSPE_robStepSplitRegSelect <- NA
      RC_robStepSplitRegSelect <- sum(unlist(robStepSplitReg_output$selections) %in% simdata$active_ind)/length(simdata$active_ind)
      PR_robStepSplitRegSelect <- sum(unlist(robStepSplitReg_output$selections) %in% simdata$active_ind)/length(unlist(robStepSplitReg_output$selections))
      
      c(MSPE_robStepSplitRegSelect, RC_robStepSplitRegSelect, PR_robStepSplitRegSelect)
      }, error = function(e){
        return(c(NA, NA, NA))
        
      }
    )
    pred_output["robStepSplitRegSelect",,i] <- robStepSplitRegSelect_final
    
    
    #robStepSplitReg with Ensemble
    
    ensemble_final <- tryCatch(
      {
        ensemble_output <- robStepSplitReg::robStepSplitReg(x = simdata$training_data$xtrain[[i]], 
                                                            y=simdata$training_data$ytrain[[i]],
                                                            model_saturation = c("fixed","p-value")[1],
                                                            n_models = n_models,
                                                            alpha = 0.05,
                                                            model_size = min(n-1, floor(p/n_models)),
                                                            compute_coef = TRUE,
                                                            enpy_opts = pense::enpy_options(retain_max = 50),
                                                            eps = 1e-1)
        MPSE_ensemble <- mean((predict(ensemble_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
        coef_ensemble_output <- coef(ensemble_output)[-1]
        RC_ensemble <- RC_PR(coef_ensemble_output, simdata$active_ind)$rc
        PR_ensemble <- RC_PR(coef_ensemble_output, simdata$active_ind)$pr
        
        c(MPSE_ensemble, RC_ensemble, PR_ensemble)
        }, error = function(e){
        return(c(NA, NA, NA))
        
      }
    )
    
    pred_output["Ensemble",,i] <- ensemble_final
    
    #robStepSplitReg Ensemble : First variables selection through step wise 
    
    ensembleSelect_final <- tryCatch(
      {
       MPSE_ensembleSelect <- NA
       RC_ensembleSelect_final <- sum(unlist(ensemble_output$selections) %in% simdata$active_ind)/length(simdata$active_ind)
       PR_ensembleSelect_final <- sum(unlist(ensemble_output$selections) %in% simdata$active_ind)/length(unlist(ensemble_output$selections))
        
       c(MPSE_ensembleSelect, RC_ensembleSelect_final, PR_ensembleSelect_final)
       }, error = function(e){
        return(c(NA, NA, NA))
      }
    )
    pred_output["EnsembleSelect",,i] <- ensembleSelect_final
    
    #HUBER
    huber_final <- tryCatch(
      {
        huber_output <- hqreg::cv.hqreg(simdata$training_data$xtrain[[i]],
                                        simdata$training_data$ytrain[[i]],
                                        alpha = 0.75, method = "huber")
        MPSE_huber <- mean((predict(huber_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
        coef_huber <- coef(huber_output)[-1]
        RC_huber <- RC_PR(coef_huber, simdata$active_ind)$rc
        PR_huber <- RC_PR(coef_huber, simdata$active_ind)$pr
        
        c(MPSE_huber, RC_huber, PR_huber)
      }, error = function(e){
        return(c(NA, NA, NA))
        
      }
    )
   
    pred_output["HuberEN",,i] <- huber_final

  }
   return (pred_output)
}

