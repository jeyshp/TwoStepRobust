#' 
#' @description generateData() produces contaminated training data and uncontaminated test data.
#' 
#' @param N Number of training sets.
#' @param n Sample size for training set.
#' @param m Sample size for test set.
#' @param p Total number of parameters.
#' @param rho Correlation within a block of active parameters.
#' @param rho.inactive Correlation between blocks of active parameters.
#' @param p.active Number of active parameters.
#' @param group.size Size of one block of active parameters.
#' @param snr Signal-to-noise ratio.
#' @param contamination.prop Contamination proportion.
#' @param contamination.scenario Casewise, cellwise marginal or cellwise correlation contamination.


generateData <- function(N, 
                         n, 
                         m, 
                         p, 
                         rho, 
                         rho.inactive, 
                         p.active, 
                         group.size, 
                         snr, 
                         contamination.prop,
                         contamination.scenario){ 
  
  # Lists to store lists
  xlist <- list()
  ylist <-list()
  
  # Casewise contamination parameters
  k_lev <- 2
  k_slo <- 100
  # Cellwise contamination parameters
  gamma <- 3
  cell.mean <- 10
  
  # Setting up correlation between and within blocks of active parameters
  sigma.mat <- matrix(0, nrow = p, ncol = p)
  sigma.mat[1:p.active, 1:p.active] <- rho.inactive
  for(group in 0:(p.active/group.size - 1))
    sigma.mat[(group*group.size+1):(group*group.size+group.size),(group*group.size+1):(group*group.size+group.size)] <- rho
  diag(sigma.mat) <- 1
  
  # Generating true coefficient vector
  trueBeta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))
  
  # Contamination of trueBeta
  beta_cont <- trueBeta
  beta_cont[trueBeta!=0] <- beta_cont[trueBeta!=0]*(1 + k_slo)
  beta_cont[trueBeta==0] <- k_slo*max(abs(trueBeta))
  
  # Noise parameter
  sigma <- as.numeric(sqrt(t(trueBeta) %*% sigma.mat %*% trueBeta)/sqrt(snr))
  
  # Simulating training data
  for(i in 1:N){

    x_train <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
    y <- x_train %*% trueBeta + rnorm(n, 0, sigma)
    
    xlist[[i]] <- x_train
    ylist[[i]] <- y

  }
  
  if(contamination.scenario == "casewise") {
    
    # Casewise Contamination 
    contamination_indices <- 1:floor(n*contamination.prop)
    for(i in 1:N) {
      for(cont_id in contamination_indices){
        
        a <- runif(p, min = -1, max = 1)
        a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
        xlist[[i]][cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
          k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
        ylist[[i]][cont_id] <- t(xlist[[i]][cont_id,]) %*% beta_cont
      }
    }
      
    } else if(contamination.scenario == "cellwise_marginal") {
      
      # Cellwise Marginal Contamination 
      for(i in 1:N) {
        
        contamination_indices <- sample(1:(n * (p + 1)), round(n * (p + 1) * contamination.prop))
        xy_train <-  cbind(xlist[[i]], ylist[[i]])
        xy_train[contamination_indices] <- NA
        for(row_id in 1:n){
          cells_id <- which(is.na(xy_train[row_id,]))
          xy_train[row_id, cells_id] <- rnorm(length(cells_id), cell.mean, 1)
        }
        xlist[[i]] <- xy_train[, -(p + 1)]
        ylist[[i]] <- xy_train[, (p + 1)]
      }
      
    } else if(contamination.scenario == "cellwise_correlation") {
      
      # Cellwise Correlation Contamination
      for(i in 1:N) {
        
        contamination_indices <- sample(1:(n * p), round(n * p * contamination.prop))
        x_train <- xlist[[i]]
        x_train[contamination_indices] <- NA
        for(row_id in 1:n){
          cells_id <- which(is.na(x_train[row_id,]))
          mu_cells <- rep(0, length(cells_id))
          sigma_cells <- sigma.mat[cells_id, cells_id]
          eigen_vec <- eigen(sigma_cells)$vectors[, length(cells_id)]
          x_train[row_id, cells_id] <- gamma * sqrt(length(cells_id)) * t(eigen_vec) /
            mahalanobis(eigen_vec, mu_cells, sigma_cells)
        }
        xlist[[i]] <- x_train
      }
      
    } else if(contamination.scenario == "mixture_marginal"){
      
      n.casewise <- floor(n*contamination.prop[1])
      
      for(i in 1:N) {
        
        # Casewise Contamination
        contamination_indices <- 1:n.casewise
        for(cont_id in contamination_indices){
          
          a <- runif(p, min = -1, max = 1)
          a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
          xlist[[i]][cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
            k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
          ylist[[i]][cont_id] <- t(xlist[[i]][cont_id,]) %*% beta_cont
        }
      
        # Cellwise Marginal Contamination 
        contamination_indices <- sample(1:((n - n.casewise) * (p + 1)), 
                                        round((n - n.casewise) * (p + 1) * contamination.prop[2]))
        xy_train <-  cbind(xlist[[i]], ylist[[i]])
        xy_train[(n.casewise + 1):n,][contamination_indices] <- NA
        for(row_id in (n.casewise + 1):n){
          cells_id <- which(is.na(xy_train[row_id,]))
          xy_train[row_id, cells_id] <- rnorm(length(cells_id), cell.mean, 1)
        }
        xlist[[i]] <- xy_train[, -(p + 1)]
        ylist[[i]] <- xy_train[, (p + 1)]
      }
    } else if(contamination.scenario == "mixture_correlation"){
      
      n.casewise <- floor(n*contamination.prop[1])
      
      for(i in 1:N) {
        
        # Casewise Contamination
        contamination_indices <- 1:n.casewise
        for(cont_id in contamination_indices){
          
          a <- runif(p, min = -1, max = 1)
          a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
          xlist[[i]][cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
            k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
          ylist[[i]][cont_id] <- t(xlist[[i]][cont_id,]) %*% beta_cont
        }
        
        # Cellwise Marginal Contamination 
        contamination_indices <- sample(1:((n - n.casewise) * p), 
                                        round((n - n.casewise) * p * contamination.prop[2]))
        x_train <-  xlist[[i]]
        x_train[(n.casewise + 1):n,][contamination_indices] <- NA
        for(row_id in (n.casewise + 1):n){
          cells_id <- which(is.na(x_train[row_id,]))
          mu_cells <- rep(0, length(cells_id))
          sigma_cells <- sigma.mat[cells_id, cells_id]
          eigen_vec <- eigen(sigma_cells)$vectors[, length(cells_id)]
          x_train[row_id, cells_id] <- gamma * sqrt(length(cells_id)) * t(eigen_vec) /
            mahalanobis(eigen_vec, mu_cells, sigma_cells)
        }
        xlist[[i]] <- x_train
      }
    }
  
  # Simulating uncontaminated test data 
  x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
  y_test <- x_test %*% trueBeta + rnorm(m, 0, sigma)
  
 return(
   list(training_data = list(xtrain = xlist, ytrain = ylist), testing_data = list(xtest = x_test, ytest = y_test), 
        pactive = p.active, n = n, sigma = sigma, active_ind = which(trueBeta != 0), p = p, trueBeta = trueBeta,
        contamination.scenario = contamination.scenario))
  
}






