#' 
#' 
#' @description generateData() produces contaminated training data and uncontaminated test data 
#' 
#' @param N Number of training sets 
#' @param n Sample size for training set
#' @param m Sample size for test set
#' @param p Total number of parameters
#' @param rho Correlation within a block of active parameters
#' @param rho.inactive Correlation between blocks of active parameters
#' @param p.active Number of active parameters
#' @param group.size Size of one block of active parameters
#' @param snr Signal to noise ratio
#' @param contamination.prop Contamination proportion


generateData <- function(N, 
                         n, 
                         m, 
                         p, 
                         rho, 
                         rho.inactive, 
                         p.active, 
                         group.size, 
                         snr, 
                         contamination.prop){ 
                         
  xlist <- list()
  ylist <-list()
  
  k_lev <- 2
  k_slo <- 100
  
  #Setting up correlation between and within blocks of active parameters
  sigma.mat <- matrix(0, nrow = p, ncol = p)
  sigma.mat[1:p.active, 1:p.active] <- rho.inactive
  for(group in 0:(p.active/group.size - 1))
    sigma.mat[(group*group.size+1):(group*group.size+group.size),(group*group.size+1):(group*group.size+group.size)] <- rho
  diag(sigma.mat) <- 1
  
  trueBeta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))
  
  #Contamination of trueBeta
  contamination_indices <- 1:floor(n*contamination.prop)
  beta_cont <- trueBeta
  beta_cont[trueBeta!=0] <- beta_cont[trueBeta!=0]*(1 + k_slo)
  beta_cont[trueBeta==0] <- k_slo*max(abs(trueBeta))
  
  sigma <- as.numeric(sqrt(t(trueBeta) %*% sigma.mat %*% trueBeta)/sqrt(snr))
  
  #Simulating and contamination of training data
  for(i in 1:N) {

    x_train <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
    y <- x_train %*% trueBeta + rnorm(n, 0, sigma)
  
    for(cont_id in contamination_indices){
      
      a <- runif(p, min = -1, max = 1)
      a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
      x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
        k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
      y[cont_id] <- t(x_train[cont_id,]) %*% beta_cont
    }
    
    xlist[[i]] <- x_train
    ylist[[i]] <- y

  }
  
 #Simulating uncontaminated test data 
  x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
  y_test <- x_test %*% trueBeta + rnorm(m, 0, sigma)
  

 return(
   list(training_data = list(xtrain = xlist, ytrain = ylist), testing_data = list(xtest = x_test, ytest = y_test), 
        pactive = p.active, n = n, sigma = sigma, active_ind = which(trueBeta != 0), p = p))
  
}






