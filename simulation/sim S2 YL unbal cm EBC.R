#==============================================================================
# Import libraries  
#==============================================================================
library(psych); library(MASS); library(mvtnorm); library(rstan); 
library(ggplot2); library(dplyr); library(foreach);
library(doParallel); library(tidyverse); library(parallel); library(mcprogress)
#==============================================================================
# Source YL codes : YL_bal_fun
#==============================================================================
source("YL_unbal_fun_cm.R")
# ==============================================================================
# Function to generate data 
#==============================================================================
# input m, A, X, beta, D
data_gen_func <- function(m, A, beta_val, D){
  # generate theta from normal(xi'beta,A)
  theta_val <- rnorm(m, mean = 0,# X%*%as.matrix(beta_val),
                     sd = sqrt(A))
  # generate y from normal(theta_i,D_i)
  y <- rnorm(m, mean = theta_val, sd = sqrt(D))
  # return theta values and y values
  data_out <- list()
  data_out$theta_val <- theta_val; data_out$y <- y
  return(data_out)
}

#==============================================================================
# Parameters: m,A,D
#==============================================================================
m_real <- 30 # m = number of small areas
A_real <- 5 # A = prior variance
# D = sampling variance
D_val <- 5*c(rep(0.01,m_real/5),rep(0.02,m_real/5),
           rep(0.1,m_real/5),rep(0.2,m_real/5),
           rep(0.3,m_real/5))
z <- 1.96

#==============================================================================
# Function for MCMC 
#==============================================================================

M_MC <- 1:100 # set the number of MCMC draws
ncores <- 5 # number of cores to be used

ST <- Sys.time()
Sys.time()
mcmc_result_EBC <- list()

mcmc_result_EBC[M_MC] <- mclapply(M_MC, function(mc_itr){
  
  #==============================================================================
  # Run MCMC
  #==============================================================================
  
  # print(paste0("MC = ",mc_itr))
  
  set.seed(mc_itr)
  
  # Step 1: Run function to generate theta and y
  data_theta_y <<- data_gen_func(m = m_real, A = A_real, 
                                 D = D_val)
  mcmc_start_time <- Sys.time()
  
  # Step 2: Compute YL function
  # print("Start YL CI")
  result_YL <- BP_func_YL(m = m_real, Di = D_val, yy = data_theta_y$y)
  
  mcmc_end_time <- Sys.time()
  
  message_parallel("MC = ",mc_itr," MCMC time: ", mcmc_end_time - mcmc_start_time)
  return(result_YL)
}, mc.cores = ncores)

Sys.time()-ST
Sys.time()


EBC <- cbind(
  # EBC
  EBC_YL = round(rowMeans(matrix(
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_EBC[[mc_itr]][, 12]))
    , byrow = F, ncol = length(M_MC))),digits = 3),
  # Length
  Length_YL = round(rowMeans(matrix(
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_EBC[[mc_itr]][, 10]))
    , byrow = F, ncol = length(M_MC))), digits = 3)
)

block_colmeans <- function(x, block_size) {
  m <- nrow(x)
  split_idx <- rep(1:ceiling(m / block_size), each = block_size)[1:m]
  
  t(sapply(split(seq_len(m), split_idx), function(idx) {
    means <- colMeans(x[idx, , drop = FALSE])
    # multiply all but last column by 100
    n <- ncol(x)
    if (n > 1) {
      means[1:(n-1)] <- 100 * means[1:(n-1)]
    }
    round(means, 1)
  }))
}
EBC_final <- block_colmeans(EBC, m_real/5)
EBC_output <- data.frame(
  YL = paste0(EBC_final[,1],"(",EBC_final[,2],")")
)


save.image("sim_bal_comb_EBC_S22_m30_1_100.RData")
write.csv(EBC_output, "sim_bal_comb_EBC_S22_m30_1_100.csv", row.names = F)

