design <- 8
# network sizes to simulate
nvec <- c(2000,3000,4000,5000,6000,10000) 
# number of simulations (not used yet)
nsim <- 1000

# other simulation parameters
TT <- 5 # time periods
K <- 4 # number of blocks
d <- 1 # dimension of latent variables
gamma <- 1.5  # persistence parameter

dmax <- 5  # largest d to try for estimation


# latent positions
latent <- matrix(cbind(.8, .3, -.6, -1.5), nrow = d, ncol = K)   

# tables to store results
table_max <- data.frame(matrix(NA, nrow = length(nvec), ncol = dmax+1))
names(table_max) <- c("n", 1:dmax)
table_mean <- data.frame(matrix(NA, nrow = length(nvec), ncol = dmax+1))
names(table_mean) <- c("n", 1:dmax)
