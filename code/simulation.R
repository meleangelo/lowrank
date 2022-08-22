# simulation
rm(list= ls())

# load functions needed for simulations
source("setup_functions.R")

# load simulation design
source("simulation_design")

asdfaf
# tables to store results
table_max <- data.frame(matrix(NA, nrow = length(nvec), ncol = length(dhatvec)+1))
names(table_max) <- c("n", paste("dhat=",dhatvec[1]), paste("dhat=",dhatvec[2]), paste("dhat=",dhatvec[3]))
table_mean <- data.frame(matrix(NA, nrow = length(nvec), ncol = length(dhatvec)+1))
names(table_mean) <- c("n", paste("dhat=",dhatvec[1]), paste("dhat=",dhatvec[2]), paste("dhat=",dhatvec[3]))



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
    Phat <- array(NA, dim = c(n,n,TT) )
    for (ii in 1:TT){
      Phat[,,ii]<- estimate_ASE(data$A[,,ii],dhat,index)
    }
    
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


