rm(list= ls())

# libraries needed
require(RSpectra)
require(irlba)
require(grdpg)
require(mclust)

# function to generate K latent positions vectors from a matrix d x K
generate_latent_posSBM <- function(latent, d, block_size){
  if (nrow(latent) != d) {
    stop("The number of rows in `latent` should equal to `d`.")
  }
  latent_positions <- matrix(rep(latent[, 1]), nrow = block_size[1], ncol = d, 
                             byrow = TRUE)
  for (k in 2:ncol(latent)) {
    latent_positions <- rbind(latent_positions, matrix(rep(latent[, k]), nrow = block_size[k], 
                                                       ncol = d, byrow = TRUE))
  }
  return(latent_positions)
}

# functoin to generate the network time series
generate_data <- function(latent_positions, gamma, time_periods, sim){
  set.seed(sim)
  n <- nrow(latent_positions)
  # generate fixed effects matrix
  alpha <- latent_positions %*% t(latent_positions)
  # generate initial matrix P_0 (stationary distribution)
  P0 <- sigmoid(alpha)/(1- sigmoid(gamma + alpha) + sigmoid(alpha))
  # generate initial adjacency A_0
  A0 <- 1*(matrix(runif(n^2), ncol = n, nrow = n) < P0 )
  A0[lower.tri(A0)] <- t(A0)[lower.tri(t(A0))]
  diag(A0) <- 0
  # initialize P and A
  P <- array(NA, dim = c(n,n,time_periods) )
  A <- array(NA, dim = c(n,n,time_periods) )
  ## observed covariates
  #X <- array(NA, dim = c(n,n,TT, nsim))
  
  #X[,,1,sim] <- matrix(rbinom(n*n , 1, .5), ncol = n, nrow = n)
  #X[,,1,sim][lower.tri(X[,,1,sim])] <- t(X[,,1,sim])[lower.tri(t(X[,,1,sim]))]
  #diag(X[,,1,sim]) <- 0
  P[,,1] <- sigmoid(gamma*A0 + alpha)
  A[,,1] <- 1*(matrix(runif(n^2), ncol = n, nrow = n) < P[,,1] )
  A[,,1][lower.tri(A[,,1])] <- t(A[,,1])[lower.tri(t(A[,,1]))]
  diag(A[,,1]) <- 0
  
  for (t in 2:time_periods){
    t1 <- t-1
    P[,,t] <- sigmoid(gamma*A[,,t1] + alpha)
    A[,,t] <- 1*(matrix(runif(n^2), ncol = n, nrow = n) < P[,,t] )
    A[,,t][lower.tri(A[,,t])] <- t(A[,,t])[lower.tri(t(A[,,t]))]
    diag(A[,,t]) <- 0
    
  }
  data <- list(A0,P0,A,P, alpha)
  names(data) <- c("A0", "P0", "A", "P", "alpha")
  return(data)
}  

# function for estimation of P with ASE period by period
estimate_ASE <- function(A_t, d, sim){
  set.seed(sim)
  # estimation of P_t
  ASE <- SpectralEmbedding(A_t, d, work  = 100)
  Ipq <- getIpq(A_t, d)
  Xhat <- ASE$X %*% sqrt(diag(ASE$D, ncol = d, nrow = d) )
  Phat <- BFcheck(Xhat %*% Ipq %*% t(Xhat))
  return(Phat)
}

# functoin to generate the network time series (NONSTATIONARY)
generate_data_nonstat <- function(latent_positions, gamma, time_periods, sim){
  set.seed(sim)
  n <- nrow(latent_positions)
  # generate fixed effects matrix
  alpha <- latent_positions %*% t(latent_positions)
  # generate initial matrix P_0 
  P0 <- sigmoid(alpha)
  # generate initial adjacency A_0
  A0 <- 1*(matrix(runif(n^2), ncol = n, nrow = n) < P0 )
  A0[lower.tri(A0)] <- t(A0)[lower.tri(t(A0))]
  diag(A0) <- 0
  # initialize P and A
  P <- array(NA, dim = c(n,n,time_periods) )
  A <- array(NA, dim = c(n,n,time_periods) )
  ## observed covariates
  #X <- array(NA, dim = c(n,n,TT, nsim))
  
  #X[,,1,sim] <- matrix(rbinom(n*n , 1, .5), ncol = n, nrow = n)
  #X[,,1,sim][lower.tri(X[,,1,sim])] <- t(X[,,1,sim])[lower.tri(t(X[,,1,sim]))]
  #diag(X[,,1,sim]) <- 0
  P[,,1] <- sigmoid(gamma*A0 + alpha)
  A[,,1] <- 1*(matrix(runif(n^2), ncol = n, nrow = n) < P[,,1] )
  diag(A[,,1]) <- 0
  
  for (t in 2:time_periods){
    t1 <- t-1
    P[,,t] <- sigmoid(gamma*A[,,t1] + alpha)
    A[,,t] <- 1*(matrix(runif(n^2), ncol = n, nrow = n) < P[,,t] )
    diag(A[,,t]) <- 0
    
  }
  data <- list(A0,P0,A,P, alpha)
  names(data) <- c("A0", "P0", "A", "P", "alpha")
  return(data)
}  

# function for estimation of X with ASE 
estimateX_ASE <- function(A_t, d, sim){
  set.seed(sim)
  # estimation of P_t
  ASE <- SpectralEmbedding(A_t, d, work  = 100)
  Xhat <- ASE$X %*% sqrt(diag(ASE$D, ncol = d, nrow = d) )
  return(Xhat)
}


# function to estimate OMNI
estimate_OMNI <- function(A, d, sim){
  set.seed(sim)
  # generate omni matrix
  n <- dim(A)[1]
  m <- dim(A)[3]
  
  Aomni <- matrix(NA, nrow = n*m , ncol = n*m )
  
  for (k in 1:m){
    for (l in 1:m) {
      Aomni[((k-1)*n+1):(k*n), ((l-1)*n+1):(l*n)] <- (A[,,k] + A[,,l])/2
    }
  }
  ASE <- SpectralEmbedding(Aomni, d, work  = 100)
  Ipq <- getIpq(Aomni, d)
  Xhat <- ASE$X %*% sqrt(diag(ASE$D, ncol = d, nrow = d) )
  Phat_omni <- Xhat %*% Ipq %*% t(Xhat)
  Phat <- array(NA, dim = c(n,n,m))
  for (k in 1:m) {
    Phat[,,k] <- Phat_omni[((k-1)*n+1):(k*n),((k-1)*n+1):(k*n)]
  }
  OMNI <- list(Xhat, Phat)
  names(OMNI) <- c("Xhat", "Phat")
  return(OMNI)
}

# function to estimate COSIE
estimate_COSIE <- function(A, d, sim){
  set.seed(sim)
  # generate omni matrix
  n <- dim(A)[1]
  m <- dim(A)[3]
  Uhat <- matrix(NA, nrow = n, ncol = d*m) # need to change this to allow multiple d's
  for (t in 1:m){
    ASE <- SpectralEmbedding(A[,,t], d, work  = 100)
    #   Ipq <- getIpq(A[,,t], d)
    Vhat_t <- ASE$X %*% sqrt(diag(ASE$D, ncol = d, nrow = d) )
    Uhat[ , ((t-1)*d+1): (t*d) ] <- Vhat_t
  }
  Vhat <- irlba(Uhat, d, work  = 100)$u
  Rhat <- array(NA, dim = c(d,d,m))
  Phat <- array(NA, dim = c(n,n,m))
  for (t in 1:m){
    Rhat[,,t] <- t(Vhat) %*% A[,,t] %*% Vhat
    Phat[,,t] <- Vhat %*% Rhat[,,t] %*% t(Vhat)
  }
  COSIE <- list(Vhat, Rhat, Phat)
  names(COSIE) <- c("Xhat", "Rhat", "Phat")
  return(COSIE)
}

estimate_GRDPG <- function(A, d, sim){
  set.seed(sim)
  n <- dim(A)[1]
  m <- dim(A)[3]
  Phat <- array(NA, dim = c(n,n,m))
  for (t in 1:m){
    Phat[,,t] <- estimate_ASE(A[,,t], d, sim)
  }
  return(Phat)
}

mean_estimate <- function(Phat,P,TT){
  n <- dim(P)[1]
  FE_est <- matrix(0, nrow = n, ncol = n)  
  for (tt in 1:TT){
    dist <- log(Phat[,,tt]/(1-Phat[,,tt] )) - log(P[,,tt]/(1-P[,,tt] ))
    FE_est <- FE_est + abs(  (1/TT)*dist  )
  }
  mean_est <- mean(FE_est)
  final <- list(FE_est,mean_est)
  names(final) <- c("FE_est","mean_est")
  return(final)
}

max_estimate <- function(Phat,P,TT){
  n <- dim(P)[1]
  FE_est <- matrix(0, nrow = n, ncol = n)  
  for (tt in 1:TT){
    dist <- log(Phat[,,tt]/(1-Phat[,,tt] )) - log(P[,,tt]/(1-P[,,tt] ))
    FE_est <- FE_est + abs(  (1/TT)*dist  )
  }
  max_est <- max(FE_est)
  final <- list(FE_est,max_est)
  names(final) <- c("FE_est","max_est")
  return(final)
}
