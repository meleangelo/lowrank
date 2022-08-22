# simulation
rm(list=ls())

# load functions needed for simulations
source("code/setup_functions.R")

# load simulation design
source("code/simulation_design.R")


# tables to store results
table_max <- data.frame(matrix(NA, nrow = length(nvec), ncol = dmax+1))
names(table_max) <- c("n", 1:dmax)
table_mean <- data.frame(matrix(NA, nrow = length(nvec), ncol = length(dhatvec)+1))
names(table_mean) <- c("n", 1:dmax)



  for (s in 1:length(nvec)){
    index = 1
    n <- nvec[s]
    cat(n, "\n")
    # Balanced case
    pi <- rep(1/K, K)
    block_size <- round(pi * n)
    blocks <- c()
    for (k in 1:length(block_size)) {
      blocks <- c(blocks, rep(k, block_size[k]))
    }
    
    # generate latent positions
    latent_pos <- generate_latent_posSBM(latent, d, block_size)
    # generate network time series data
    #data <- generate_data(latent_pos, gamma, TT, index)
    data <- generate_data_nonstat(latent_pos, gamma, TT, index)
    
    # estimate Phat
    Xhat <- array(NA, dim = c(n,dmax,TT) )
    for (tt in 1:TT){
      Xhat[,,tt]<- estimateX_ASE(data$A[,,tt],dmax,index)
    }
    
    # add Ipq and dim select here !!!
    
    # compute max and mean distance 
    diff_max <- max_estimate(Phat,data$P,TT)
    diff_mean <- mean_estimate(Phat,data$P,TT)
    
    # record in tables
    table_max[,1] <- nvec
    table_max[s,dd+1] <- diff_max$max_est
    table_mean[,1] <- nvec
    table_mean[s,dd+1] <- diff_mean$mean_est
  }

table_max
table_mean


