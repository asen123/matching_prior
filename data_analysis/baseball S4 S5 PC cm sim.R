# Using common mean with 
  # 1. arcsine transformed simulated y
  # 2. simulated normal y

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

#==============================================================================
# Real data analysis : Baseball data
#==============================================================================
# Columns:
# "player": name of players;
# "est_dir": batting average of player d for first 45 at-bats;
# "x1": batting average of player d through the end of 1969 season;
# "x2": total number of times player d at-bats through the end of 1969 season;
# "actual": true season batting average of player d for the 1970 season;
baseball = read.table("baseball.txt")
colnames(baseball) = c("player", "est_dir", "x1", "x2", "actual")
m_real <- as.numeric(nrow(baseball))

# arc sine transformation
arcsin <- function(x) return(sqrt(45)*asin(2*x-1))


## Simulate binomial from the actual and arcsine transformation (S4)
set.seed(1234)
baseball$y <- arcsin(rbinom(m_real, 45, baseball$actual)/45)
## Theta value
baseball$theta <- arcsin(baseball$actual)
## Simulate normal from with mean theta (S5)
# set.seed(1234)
# baseball$y <- rnorm(m_real, baseball$theta, 1)

# common mean
mu <- -3.275
z <- 1.96 
D_val <- rep(1,m_real) 
A_real <- 0.26582 # A = prior variance

data_theta_y <- list(theta = baseball$theta -mu , 
                     y = baseball$y - mu)

D_val <- rep(1,m_real) 

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
M_MC <- 1:100
post_cov = round(100*cbind(
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
),1)



save.image("baseball_pb_S5_1_100.RData")
write.csv(post_cov, "baseball_pb_S5_1_100.csv", row.names = F)

# save.image("baseball_pb_S6_1_100.RData")
# write.csv(post_cov, "baseball_pb_S6_1_100.csv", row.names = F)

