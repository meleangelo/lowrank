for (dd in 1:dmax){
  dhat <- dd
  for (s in 1:length(nvec)){
    index = 1
    n <- nvec[s]
    cat("Design", design, "; n=",n, ", dhat=", dhat, "Started at:", Sys.time() , "\n")
    # Balanced blocks
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
    table_max[s,1] <- nvec[s]
    table_max[s,dd+1] <- diff_max$max_est
    table_mean[s,1] <- nvec[s]
    table_mean[s,dd+1] <- diff_mean$mean_est
  }
}

# save file
filename <- paste("output/simulation_design", design, ".RData", sep = "")
save( table_max, table_mean, file = filename)

# tex tables
table_filename <- paste("output/tables_design", design, ".tex", sep = "")
tex_max <- xtable(table_max, digits = c(0,0,rep(4,dmax)))
tex_mean <- xtable(table_mean, digits = c(0,0,rep(4,dmax)))
print.xtable(tex_max, file = table_filename, include.rownames = FALSE )
print.xtable(tex_mean, file = table_filename, include.rownames = FALSE , append = TRUE )
