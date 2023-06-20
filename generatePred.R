
#generatePred() fits a model using PENSE and sparseLTS and stores the MSPE,specificity and sensitivity of each training set for each model

source("RC_PR.R")
generatePred <- function(simdata,n_models, ...) {
  
  N = length(simdata$training_data$xtrain) #size of training dataset
  p.active = simdata$pactive # number of active predictors
  n = simdata$n
  
  #test data
  xtestdata = simdata$testing_data$xtest 
  ytestdata = simdata$testing_data$ytest
  
  
  pred_output <- array(dim = c(6,3,N))
  colnames(pred_output) = c("MSPE","RC", "PR")
  rownames(pred_output) = c("EN", "Pense", "sparseLTS", "robStepSplitReg","Ensemble", "HuberEN")
  
  
  
  for(i in 1:N) {
    
    en_output <- glmnet::cv.glmnet(x = simdata$training_data$xtrain[[i]],
                                   y = simdata$training_data$ytrain[[i]],
                                   alpha = 0.75)
    pred_output["EN", "MSPE", i] <- mean((predict(en_output,xtestdata, lambda = "lambda.min")- ytestdata)^2)/simdata$sigma^2
    coef_en <- coef(en_output,lambda = "lambda.min")[-1]
    pred_output["EN", "RC", i] <- RC_PR(coef_en, simdata$active_ind)$rc
    pred_output["EN", "PR", i] <- RC_PR(coef_en, simdata$active_ind)$pr
    
    pense_output <- pense::adapense_cv(x = simdata$training_data$xtrain[[i]],
                                       y=simdata$training_data$ytrain[[i]],
                                       alpha = 0.75, cv_k = 10, cv_repl = 1, ...)
    pred_output["Pense", "MSPE", i] <- mean((predict(pense_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
    coef_pense <- coef(pense_output)[-1]
    pred_output["Pense", "RC", i] <- RC_PR(coef_pense, simdata$active_ind)$rc
    pred_output["Pense", "PR", i] <- RC_PR(coef_pense, simdata$active_ind)$pr
    
    lambda_grid <- lambda0(simdata$training_data$xtrain[[i]], simdata$training_data$ytrain[[i]])
    sparseLTS_output <- robustHD::sparseLTS(x = simdata$training_data$xtrain[[i]], 
                                            y = simdata$training_data$ytrain[[i]], 
                                            lambda =lambda_grid)
    pred_output["sparseLTS", "MSPE" , i] <- mean((predict(sparseLTS_output, xtestdata) - ytestdata)^2)/simdata$sigma^2
    coef_sparseLTS <- coef(sparseLTS_output)[-1]
    pred_output["sparseLTS", "RC", i] <- RC_PR(coef_sparseLTS, simdata$active_ind)$rc
    pred_output["sparseLTS", "PR", i] <- RC_PR(coef_sparseLTS, simdata$active_ind)$pr
    
    robStepSplitReg_output <- robStepSplitReg::robStepSplitReg(x = simdata$training_data$xtrain[[i]], 
                                                               y=simdata$training_data$ytrain[[i]],
                                                               model_saturation = c("fixed","p-value")[1],
                                                               alpha = 0.05,
                                                               model_size = n-1,
                                                               compute_coef = TRUE)
    pred_output["robStepSplitReg", "MSPE", i] <- mean((predict(robStepSplitReg_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
    coef_robStepSplitReg <- coef(robStepSplitReg_output)[-1]
    pred_output["robStepSplitReg", "RC", i] <- RC_PR(coef_robStepSplitReg, simdata$active_ind)$rc
    pred_output["robStepSplitReg", "PR", i] <- RC_PR(coef_robStepSplitReg, simdata$active_ind)$pr
    
    ensemble_output <- robStepSplitReg::robStepSplitReg(x = simdata$training_data$xtrain[[i]], 
                                                               y=simdata$training_data$ytrain[[i]],
                                                               model_saturation = c("fixed","p-value")[1],
                                                               n_models = n_models,
                                                               alpha = 0.05,
                                                               model_size = n-1,
                                                               compute_coef = TRUE)
    pred_output["Ensemble", "MSPE", i] <- mean((predict(ensemble_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
    coef_ensemble_output <- coef(ensemble_output)[-1]
    pred_output["Ensemble", "RC", i] <- RC_PR(coef_ensemble_output, simdata$active_ind)$rc
    pred_output["Ensemble", "PR", i] <- RC_PR(coef_ensemble_output, simdata$active_ind)$pr
    
    huber_output <- hqreg::cv.hqreg(simdata$training_data$xtrain[[i]],
                                    simdata$training_data$ytrain[[i]],
                                    alpha = 0.75, method = "huber")
    pred_output["HuberEN", "MSPE", i] <- mean((predict(huber_output,xtestdata)- ytestdata)^2)/simdata$sigma^2
    coef_huber <- coef(huber_output)[-1]
    pred_output["HuberEN", "RC", i] <- RC_PR(coef_huber, simdata$active_ind)$rc
    pred_output["HuberEN", "PR", i] <- RC_PR(coef_huber, simdata$active_ind)$pr
    

  }
   return (pred_output)
}

