#==============================================================================
# Import libraries  
#==============================================================================
library(psych); library(MASS); library(mvtnorm); library(rstan); 
library(ggplot2); library(dplyr); library(foreach);
library(doParallel); library(tidyverse); library(parallel); library(mcprogress)

# ==============================================================================
source("AS_bal_fun.R")
source("YL_bal_fun.R")
# ==============================================================================
# Function to generate data 
#==============================================================================
# input m, A, X, beta, D
data_gen_func <- function(m, A, X, beta_val, D){
  # generate theta from normal(xi'beta,A)
  theta_val <- rnorm(m, mean = X%*%as.matrix(beta_val), 
                     sd = sqrt(A))
  # generate y from normal(theta_i,D_i)
  y <- rnorm(m, mean = theta_val, sd = sqrt(D))
  # return theta values and y values
  data_out <- list()
  data_out$theta_val <- theta_val; data_out$y <- y
  return(data_out)
}
#==============================================================================
# Parameters: m,A,D,X,beta 
#==============================================================================
m_real <- 30 # m = number of small areas
A_real <- 10 # A = prior variance
# D = sampling variance
D_val <- rep(5,m_real)
# beta vector (beta0,beta1,beta2)
p <- 3 # dimension
beta_real <- c(2,-5,8)
# X matrix with 1st column as 1 for beta0 and 2 other columns 
# for beta1, beta2 (change accordingly)
set.seed(5675)
X_mat <- cbind(rep(1,m_real), 
               matrix(rnorm(m_real*(p-1),mean=0,sd=1),
                      nrow = m_real, ncol=(p-1)))
z <- 1.96
p_val <- ncol(X_mat)

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
  data_theta_y <<- data_gen_func(m = m_real, A = A_real, D = D_val,
                                 X = X_mat, beta_val = beta_real)
  
  mcmc_start_time <- Sys.time()
  
  # Step 2: Compute YL function
  # print("Start YL CI")
  result_YL <- EBLUP_func_AS(m = m_real, Di = D_val, 
                             X = X_mat, yy = data_theta_y$y)
  result_AS <- EBLUP_func_YL(m = m_real, Di = D_val, 
                             X = X_mat, yy = data_theta_y$y)
  
  mcmc_end_time <- Sys.time()
  
  message_parallel("MC = ",mc_itr," MCMC time: ", mcmc_end_time - mcmc_start_time)
  return(list(result_YL,result_AS))
}, mc.cores = ncores)

Sys.time()-ST
Sys.time()


EBC <- cbind(
  # EBC
  EBC_AS = round(rowMeans(matrix(
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_EBC[[mc_itr]][[1]][, 12]))
    , byrow = F, ncol = length(M_MC))),digits = 3),
  EBC_YL = round(rowMeans(matrix(
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_EBC[[mc_itr]][[2]][, 12]))
    , byrow = F, ncol = length(M_MC))),digits = 3),
  # Length
  Length_AS = round(rowMeans(matrix(
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_EBC[[mc_itr]][[1]][, 10]))
    , byrow = F, ncol = length(M_MC))), digits = 3),
  Length_YL = round(rowMeans(matrix(
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_EBC[[mc_itr]][[2]][, 10]))
    , byrow = F, ncol = length(M_MC))), digits = 3)
)

block_colmeans <- function(x, block_size) {
  m <- nrow(x)
  split_idx <- rep(1:ceiling(m / block_size), each = block_size)[1:m]
  
  t(sapply(split(seq_len(m), split_idx), function(idx) {
    means <- colMeans(x[idx, , drop = FALSE])
    # multiply all but last two columns by 100
    n <- ncol(x)
    if (n > 2) {
      means[1:(n-2)] <- 100 * means[1:(n-2)]
    }
    round(means, 1)
  }))
}
EBC_final <- block_colmeans(EBC, m_real/5)
EBC_output <- data.frame(AS = paste0(EBC_final[,1],"(",EBC_final[,3],")"),
                         YL = paste0(EBC_final[,2],"(",EBC_final[,4],")")
)


save.image("sim_bal_comb_EBC_S33_m30_1_100.RData")
write.csv(EBC_output, "sim_bal_comb_EBC_S33_m30_1_100.csv", row.names = F)

