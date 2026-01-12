#==============================================================================
# Import libraries  
#==============================================================================
library(psych); library(MASS); library(mvtnorm); library(rstan); 
library(ggplot2); library(dplyr); library(foreach);
library(doParallel); library(tidyverse); library(parallel); library(mcprogress)
#==============================================================================
# Source YL codes : YL_unbal_fun_cm
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
# Parameters: m,A,D,X,beta 
#==============================================================================
m_real <- 30 # m = number of small areas
A_real <- 1 # A = prior variance
# D = sampling variance
D_val <- c(rep(0.01,m_real/5),rep(0.02,m_real/5),
           rep(0.1,m_real/5),rep(0.2,m_real/5),
           rep(0.3,m_real/5))
z <- 1.96
data_theta_y <- data_gen_func(m = m_real, A = A_real, 
                              D = D_val)

# Run function to generate I_YL CI
result_YL <- BP_func_YL(m = m_real, Di = D_val, yy = data_theta_y$y)

#==============================================================================
# Function for MCMC for AS prior
#==============================================================================

mcmc_result_comb <- list()

M_MC <- 51:100 # set the number of MCMC draws
ncores <- 5 # number of cores to be used

ST <- Sys.time()
Sys.time()

mcmc_result_comb[M_MC] <- mclapply(M_MC, function(mc_itr){
  
  #==============================================================================
  # Define stan programs for different priors: data, parameters and model
  #==============================================================================
  
  stan.model.custom.AS  <- "
  data {
    int<lower=1> m;
    vector[m] y;
    real<lower=0> D[m];
    int<lower=1,upper=m> i; // index for pi_i(A)
  }

  parameters {
    real<lower=0> A;
    vector[m] theta;
  }

  model {
    theta ~ normal(0, sqrt(A));  // Prior for theta_i | beta, A
    y ~ normal(theta, sqrt(D));   // Likelihood: y_i | theta_i
  
    // Custom prior for A
  {
      vector[m] denom;
      for(j in 1:m)
        denom[j] = A + D[j];
    
      target += log(sum(inv(square(denom))));       //  sum 1 / (A + D_u)^2
      target += 2 * log(A + D[i]);             // (A + D_i)^2
      target += log(A);                        // A
      target += (y[i])^2/(A + D[i]);
    }
  }
  "
  
  #==============================================================================
  func_AS <- function(data_theta_y, i_val){
    # print(paste0("start of area ",i_val))
    
    # options(buildtools.check=function(action) TRUE)
    # input data for area dependent prior
    data.stan <- list(m = m_real, i = i_val, 
                      y = data_theta_y$y, D = D_val)
    
    prob.stanfit = stan(model_code = stan.model.custom.AS, 
                        data = data.stan,
                        chains = 4, iter = 1000, 
                        warmup = 500, refresh = 0)
    
    summ.stanfit <- summary(prob.stanfit, pars = c("theta"), 
                            probs = c(0.025, 0.975), digits_summary = 4)
    MC.summary <- data.frame(
      post_mean = summ.stanfit$summary[i_val,1], # posterior mean
      post_var = (summ.stanfit$summary[i_val,3])^2, # posterior variance
      lower = summ.stanfit$summary[i_val,4], # lower limit of CI
      upper = summ.stanfit$summary[i_val,5]) # upper limit of CI
    
    # summary of A values
    ex.post = rstan::extract(prob.stanfit)
    # print(summary(ex.post$A))
    # posterior coverage: if theta sample is inside YL and AS
    MC.summary$YL_lower <-result_YL$YL_lower[i_val]
    MC.summary$YL_upper <- result_YL$YL_upper[i_val]
    
    # MC.summary$AS_lower <-result_AS$AS_lower[i_val]
    # MC.summary$AS_upper <- result_AS$AS_upper[i_val]
    
    theta_sample <- as.vector(ex.post$theta[,i_val])
    
    MC.summary$chk_YL <- mean(ifelse(theta_sample >= MC.summary$YL_lower & 
                                       theta_sample <= MC.summary$YL_upper, 1,0))
    # MC.summary$chk_AS <- mean(ifelse(theta_sample >= MC.summary$AS_lower & 
    #                                    theta_sample <= MC.summary$AS_upper, 1,0))
    return(list(MC.summary, theta_sample))
  }
  
  #==============================================================================
  # Run MCMC
  
  # print(paste0("MC = ",mc_itr))
  
  set.seed(mc_itr)
  
  mcmc_start_time <- Sys.time()
  
  # Run stan programs
  result_final <- lapply(1:m_real, function(i_val) func_AS(data_theta_y, i_val))
  
  mcmc_end_time <- Sys.time()
  
  message_parallel("MC = ",mc_itr," MCMC time: ", mcmc_end_time - mcmc_start_time)
  return(result_final)
}, mc.cores = ncores)

Sys.time()-ST
Sys.time()


# ==============================================================================

post_cov = cbind(
  # AS_cov = rowMeans(matrix(unlist(lapply(1:m_real, function(i){
  #   unlist(lapply(M_MC, function(mc_itr) mcmc_result_comb[[mc_itr]][[i]][[1]]$chk_AS))
  # } )), byrow = T, ncol = length(M_MC))),
  YL_cov = rowMeans(matrix(unlist(lapply(1:m_real, function(i){
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_comb[[mc_itr]][[i]][[1]]$chk_YL))
  } )), byrow = T, ncol = length(M_MC))),
  # med_AS = unlist(apply(matrix(unlist(lapply(1:m_real, function(i){
  #   unlist(lapply(M_MC, function(mc_itr) mcmc_result_comb[[mc_itr]][[i]][[1]]$chk_AS))
  # } )), byrow = T, ncol = length(M_MC)),1,median)),
  med_YL = unlist(apply(matrix(unlist(lapply(1:m_real, function(i){
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_comb[[mc_itr]][[i]][[1]]$chk_YL))
  } )), byrow = T, ncol = length(M_MC)),1,median))
)



block_colmeans <- function(x, block_size) {
  m <- nrow(x)
  split_idx <- rep(1:ceiling(m / block_size), each = block_size)[1:m]
  t(sapply(split(seq_len(m), split_idx), function(idx) round(100*colMeans(x[idx, , drop = FALSE]),1)))
}

post_cov <- block_colmeans(post_cov, m_real/5)

save.image("sim_unbal_cm_post_cov_S22_m30_1_100.RData")
write.csv(post_cov, "sim_unbal_cm_post_cov_S22_m30_1_10 0.csv", row.names = F)
