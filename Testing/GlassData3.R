
#Required packages and source files 
install.packages("ltsspca")
library(ltsspca)
data("glass")

source("SparseShootingS/sparseShootingS.R")

#Checking for NA values 
Glass_df <- data.frame(Glass)
NA_val <- sum(is.na(Glass_df))


# Initialization of output array
response <- 4 #Number of response variables 
n_models <- 5 #Number of models for ensemble
p <- 737 #Number of predictors
N <- 50
n <- 30 #training sample size

pred_output <- array(dim = c(10, 2, N, response))
colnames(pred_output) <- c("MSPE", "CPU")
rownames(pred_output) <- c("ElasticNet", "DDC_ElasticNet", 
                           "sparseShooting",
                           "Pense", "sparseLTS",
                           "HuberEN",
                           "robStepSplitReg", "robStepSplitRegEnsemble", 
                           "srlars","srlarsEnsemble")
dimnames(pred_output)[[3]] <- c(seq(1:N))
dimnames(pred_output)[[4]] <- c("Na2O","MgO","Al2O3","SiO2","P2O5","SO3","Cl","K2O","CaO","MnO","Fe2O3","BaO","PbO")[c(1,4,8,9)]

response_output <- array(dim = c(10, 2, N))
colnames(response_output) <- c("MSPE", "CPU")
rownames(response_output) <- c("ElasticNet", "DDC_ElasticNet", 
                           "sparseShooting",
                           "Pense", "sparseLTS",
                           "HuberEN",
                           "robStepSplitReg", "robStepSplitRegEnsemble", 
                           "srlars","srlarsEnsemble")

#Removing first 13 predictors and shuffling data 
full_data <- list(x = Glass_df[,14:750], y = y[,c(1,4,8,9)])
#full_data$y = scale(full_data$y)
shuffle_id <- sample(1:nrow(full_data$x))
full_data$x[shuffle_id,]
full_data$y[shuffle_id,]

#Choosing 30 training ids

train_ids <- lapply(1:N, function(x)
  return(sample(1:180, 30)))

filename_total <- paste0("results_for_total_responses_N=",N,"_.RData")

for (r in 1:response) {
  # Outliers in response variables 
  outliers <- cellWise::DDC(cbind(full_data$x, full_data$y[,r]),
                            DDCpars = list(fastDDC = TRUE, silent = TRUE))
  outliers_mat <- matrix(0, nrow = 180, ncol = 738)
  outliers_mat[outliers$indcells] <- 1
  
  y_outliers_mat <- outliers_mat[,738]
  outliers_id <- which(y_outliers_mat == 1)
  
  filename <- paste0("results_for_",dimnames(pred_output)[[4]][[r]],"_N=",N,"_.RData")
  
  for (i in 1:N){
  #Splitting data into training and testing data set
    #Training Data
    xtraindata <- full_data$x[train_ids[[i]],]
    ytraindata <- full_data$y[train_ids[[i]],r]
    
    # Test data with contaminated samples removed
    test_train_ids <- c(train_ids[[i]], outliers_id)
    test_train_ids <- test_train_ids[!duplicated(test_train_ids)]
    
    xtestdata <- full_data$x[-test_train_ids,]
    ytestdata <- full_data$y[-test_train_ids,r]

    
    # Elastic net
    en_final <- tryCatch({
      en_cpu <- system.time(en_output <-   glmnet::cv.glmnet(x = as.matrix(xtraindata),
                                                             y = ytraindata,
                                                             alpha = 3/4))
      
      MSPE_en <- mean((predict(en_output,as.matrix(xtestdata))- ytestdata)^2)
      CPU_en <- en_cpu["elapsed"]
      
      c(MSPE_en, CPU_en)
      
    }, error = function(e){
      return(c(NA, NA))
    })
    
    pred_output["ElasticNet",,i,r] <- en_final
    response_output["ElasticNet",,i] <- en_final

    # sparseshootingS_final <- tryCatch({
    # 
    #   sparseS_cpu <- system.time(sparseS_output <- sparseshooting(x =  as.matrix(xtraindata),
    #                                                               y =  ytraindata,
    #                                                               wvalue = 3, nlambda = 100))
    # 
    #   sparseS_preds <- sparseS_output$coef[1] + as.matrix(xtestdata) %*% sparseS_output$coef[-1]
    #   MSPE_sparseS <- mean((sparseS_preds -ytestdata)^2)
    #   CPU_sparseS <- sparseS_cpu["elapsed"]
    # 
    #   c(MSPE_sparseS, CPU_sparseS)
    # 
    # }, error = function(e){
    #   return(c(NA, NA))
    # })
    # pred_output["sparseShooting",, i,r] <- sparseshootingS_final
    # response_output["sparseShooting",, i] <- sparseshootingS_final


    DDC_ElasticNet_final <- tryCatch({
      
      DDCxy <- cellWise::DDC(cbind(xtraindata,ytraindata),
                             DDCpars = list(fastDDC = TRUE, silent = TRUE))
      x_imp <- DDCxy$Ximp[, -ncol(DDCxy$Ximp)]
      y_imp <- DDCxy$Ximp[, ncol(DDCxy$Ximp)]
      
      DDC_en_cpu <- system.time(enDDC_output <- glmnet::cv.glmnet(x = x_imp,
                                                                  y = y_imp,
                                                                  alpha = 0.75))
      MSPE_DDC_ElasticNet<- mean((predict(enDDC_output, as.matrix(xtestdata), lambda = "lambda.min")- ytestdata)^2)
      
      
      CPU_DDC_ElasticNet <- DDC_en_cpu["elapsed"]
      
      c(MSPE_DDC_ElasticNet,CPU_DDC_ElasticNet)
      
    }, error = function(e){
      return(c(NA, NA))
    })
    pred_output["DDC_ElasticNet",,i,r] <- DDC_ElasticNet_final
    response_output["DDC_ElasticNet",,i] <- DDC_ElasticNet_final
    
    
    #PENSE
    
    # pense_final <- tryCatch({
    # 
    #   pense_cpu <- system.time(pense_output <- pense::adapense_cv(x = xtraindata[,-c(1:14)],
    #                                                               y= ytraindata,
    #                                                               alpha = 0.75, cv_k = 10, cv_repl = 1))
    #   MSPE_pense <- mean(unlist(predict(pense_output,xtestdata[,-c(1:14)])- ytestdata)^2)
    #   CPU_pense <- pense_cpu["elapsed"]
    # 
    #   c(MSPE_pense)
    # }, error = function(e){
    #   return(c(NA, NA))
    # })
    # pred_output["Pense",i, 1] <- pense_final
    
    
    
    # sparseLTS_final <- tryCatch({
    # 
    #   lambda_max <- robustHD::lambda0(as.matrix(xtraindata), ytraindata)
    #   lambda_grid = rev(exp(seq(log(1e-2*lambda_max),log(lambda_max),length = 50)))
    #   sparseLTS_cpu <- system.time(sparseLTS_output <- robustHD::sparseLTS(x = as.matrix(xtraindata),
    #                                                                        y = ytraindata,
    #                                                                        lambda = lambda_grid,
    #                                                                        mode = "lambda",
    #                                                                        tol = 1e-2))
    #   MSPE_sparseLTS <- mean((predict(sparseLTS_output, as.matrix(xtestdata)) - ytestdata)^2)
    # 
    #   CPU_sparseLTS <- sparseLTS_cpu["elapsed"]
    # 
    #   c(MSPE_sparseLTS, CPU_sparseLTS)
    # 
    # }, error = function(e){
    #   return(c(NA, NA))
    # })
    # pred_output["sparseLTS",, i,r] <- sparseLTS_final
    # response_output["sparseLTS",, i] <- sparseLTS_final

    
    # robStepSplitReg
    robStepSplitReg_final <- tryCatch({
      
      robStepSplitReg_cpu <- system.time(robStepSplitReg_output <- robStepSplitReg::robStepSplitReg(x =as.matrix(xtraindata), 
                                                                                                    y= ytraindata,
                                                                                                    model_saturation = c("fixed","p-value")[1],
                                                                                                    alpha = 0.05,
                                                                                                    model_size = n-1,
                                                                                                    compute_coef = TRUE))
      
      MPSE_robStepSplitReg <- mean((predict(robStepSplitReg_output,as.matrix(xtestdata))- ytestdata)^2)
      CPU_robStepSplitReg <- robStepSplitReg_cpu["elapsed"]
      
      c(MPSE_robStepSplitReg,CPU_robStepSplitReg)
      
    }, error = function(e){
      return(c(NA, NA))
      
    })
    pred_output["robStepSplitReg",,i,r] <- robStepSplitReg_final
    response_output["robStepSplitReg",,i] <- robStepSplitReg_final
    
    
    # robStepSplitReg with Ensemble
    
    robStepSplitReg_ensemble_final <- tryCatch({
      
      robStepSplitReg_ensemble_cpu <- system.time(robStepSplitReg_ensemble_output <- robStepSplitReg::robStepSplitReg(x = as.matrix(xtraindata), 
                                                                                                                      y= ytraindata,
                                                                                                                      model_saturation = c("fixed","p-value")[1],
                                                                                                                      n_models = n_models,
                                                                                                                      alpha = 0.05,
                                                                                                                      model_size = min(n-1, floor(p/n_models)),
                                                                                                                      compute_coef = TRUE))
      MPSE_robStepSplitReg_ensemble <- mean((predict(robStepSplitReg_ensemble_output,as.matrix(xtestdata))- ytestdata)^2)
      CPU_robStepSplitReg_ensemble <- robStepSplitReg_ensemble_cpu["elapsed"]
      
      c(MPSE_robStepSplitReg_ensemble,CPU_robStepSplitReg_ensemble)
      
    }, error = function(e){
      return(c(NA, NA))
    })
    pred_output["robStepSplitRegEnsemble",,i,r] <- robStepSplitReg_ensemble_final
    response_output["robStepSplitRegEnsemble",,i] <- robStepSplitReg_ensemble_final
    
    # Split Robust LARS
    srlars_final <- tryCatch({
      
      srlars_cpu <- system.time(srlars_output <- srlars::srlars(as.matrix(xtraindata),
                                                                ytraindata, 
                                                                n_models = 1,
                                                                model_saturation = c("fixed", "p-value")[1],
                                                                alpha = 0.05, model_size = n-1,
                                                                robust = TRUE,
                                                                compute_coef = TRUE,
                                                                en_alpha = 1/4))
      
      MSPE_srlars<- mean((predict(srlars_output,as.matrix(xtestdata))- ytestdata)^2)
      CPU_srlars  <- srlars_cpu["elapsed"]
      
      c(MSPE_srlars,CPU_srlars)
      
    }, error = function(e){
      return(c(NA, NA))
    })
    pred_output["srlars",,i,r] <- srlars_final
    response_output["srlars",,i] <- srlars_final
    
    
    # Split Robust LARS Ensemble
    srlarsEnsemble_final <- tryCatch({
      
      
      srlarsEnsemble_cpu <- system.time(srlars_outputEnsemble <- srlars::srlars(as.matrix(xtraindata),
                                                                                ytraindata, 
                                                                                n_models = n_models,
                                                                                model_saturation = c("fixed", "p-value")[1],
                                                                                alpha = 0.05, model_size = min(n-1, floor(p/n_models)),
                                                                                robust = TRUE,
                                                                                compute_coef = TRUE,
                                                                                en_alpha = 1/4))
      
      MSPE_srlarsEnsemble<- mean((predict(srlars_outputEnsemble,as.matrix(xtestdata))- ytestdata)^2)
      CPU_srlarsEnsemble  <- srlarsEnsemble_cpu["elapsed"]
      
      c(MSPE_srlarsEnsemble,CPU_srlarsEnsemble)
      
    }, error = function(e){
      return(c(NA, NA))
    })
    pred_output["srlarsEnsemble",,i,r] <- srlarsEnsemble_final
    response_output["srlarsEnsemble",,i] <- srlarsEnsemble_final
    
    
    # Huber-EN
    huber_final <- tryCatch({
      
      HUBER_cpu <- system.time(huber_output <- hqreg::cv.hqreg(as.matrix(xtraindata),
                                                               ytraindata,
                                                               alpha = 0.75, method = "huber"))
      MPSE_huber <- mean((predict(huber_output,as.matrix(xtestdata))- ytestdata)^2)
      CPU_huber  <- HUBER_cpu["elapsed"]
      
      c(MPSE_huber,CPU_huber)
      
    }, error = function(e){
      return(c(NA, NA))
    })
    pred_output["HuberEN",,i,r] <- huber_final
    response_output["HuberEN",,i] <- huber_final
    }
  #save(response_output, file = filename)
}

save(pred_output, file = filename_total)


compounds <- apply(pred_output, c(1,2,4),mean, na.rm = TRUE)
ranks_scaled <- sapply(1:4, function(x, compounds) return (rank(compounds[,,x][,1])) ,compounds = compounds)
keep <- which(apply(y, 2, var) >1)
apply(ranks_scaled,1, mean)

compounds <- apply(pred_output, c(1,2,4),mean)
ranks_unscaled <- sapply(1:13, function(x, compounds) return (rank(compounds[,,x][,1])) ,compounds = compounds)
apply(ranks_unscaled[,keep],1, mean)




