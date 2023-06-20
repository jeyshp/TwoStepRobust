
#generateData() produces a training data set and a test data set

generateData <- function(N, #number of training sets
                         n, #sample size for training set
                         m, #sample size for test set
                         p, #number of peters
                         rho, #correlation within a block of peters
                         rho.inactive, #correlation between blocks 
                         p.active, #number of active peters
                         group.size, #size of one block
                         snr, 
                         contamination.prop){  #contamination proportion
                         
  xlist <- list()
  ylist <-list()
  
  k_lev <- 2
  k_slo <- 100
  
  sigma.mat <- matrix(0, nrow = p, ncol = p)
  sigma.mat[1:p.active, 1:p.active] <- rho.inactive
  for(group in 0:(p.active/group.size - 1))
    sigma.mat[(group*group.size+1):(group*group.size+group.size),(group*group.size+1):(group*group.size+group.size)] <- rho
  diag(sigma.mat) <- 1
  
  trueBeta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))
  
  #contamination of trueBeta
  contamination_indices <- 1:floor(n*contamination.prop)
  beta_cont <- trueBeta
  beta_cont[trueBeta!=0] <- beta_cont[trueBeta!=0]*(1 + k_slo)
  beta_cont[trueBeta==0] <- k_slo*max(abs(trueBeta))
  
  sigma <- as.numeric(sqrt(t(trueBeta) %*% sigma.mat %*% trueBeta)/sqrt(snr))
  
  #training data + contamination
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
  
 #test data 
  x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
  y_test <- x_test %*% trueBeta + rnorm(m, 0, sigma)
  

 return(
   list(training_data = list(xtrain = xlist, ytrain = ylist), testing_data = list(xtest = x_test, ytest = y_test), 
        pactive = p.active, n = n, sigma = sigma, active_ind = which(trueBeta != 0)))
  
}






