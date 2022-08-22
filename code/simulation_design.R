
# network sizes to simulate
nvec <- c(2000,3000,4000,5000,6000,10000) 
# number of simulations (not used yet)
nsim <- 1000

# other simulation parameters
TT <- 5 # time periods
K <- 2 # number of blocks
d <- 1 # dimension of latent variables
gamma <- .1  # persistence parameter

dmax <- 5  # largest d to try for estimation


# latent positions
latent <- matrix(cbind(.8, -1.5), nrow = d, ncol = K)   

