#==============================================================================
# Import libraries  
#==============================================================================
library(psych); library(MASS); library(mvtnorm); library(rstan); 
library(ggplot2); library(dplyr); library(foreach);
library(doParallel); library(tidyverse); library(parallel); library(mcprogress)
#==============================================================================
# Source YL code : YL_unbal_fun, AS_unbal_fun
#==============================================================================
source("AS_unbal_fun.R")
source("YL_unbal_fun.R")
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
m_real <- 15 # m = number of small areas
A_real <- 1*5 # A = prior variance
# D = sampling variance
D_val <- c(rep(0.01*5,m_real/5),rep(0.02*5,m_real/5),
           rep(0.1*5,m_real/5),rep(0.2*5,m_real/5),
           rep(0.3*5,m_real/5))
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
data_theta_y <- data_gen_func(m = m_real, A = A_real, 
                              X = X_mat, beta_val = beta_real, 
                              D = D_val)

# result_flat <- result_YL <- result_AS <- list()
# ST <- Sys.time()
# Run function to generate I_YL CI
result_YL <- EBLUP_func_YL(m = m_real, Di = D_val, 
                        X = X_mat, yy = data_theta_y$y)
result_AS <- EBLUP_func_AS(m = m_real, Di = D_val, 
                           X = X_mat, yy = data_theta_y$y)

#==============================================================================
# Function for MCMC for AS prior
#==============================================================================

M_MC <- 1:100 # set the number of MCMC draws
ncores <- 5 # number of cores to be used

ST <- Sys.time()
Sys.time()
mcmc_result_comb <- list()

mcmc_result_comb[M_MC] <- mclapply(M_MC, function(mc_itr){
  
  #==============================================================================
  # Define stan programs for different priors: data, parameters and model
  #==============================================================================
  
  stan.model.custom.AS  <- "
  data {
    int<lower=1> m;
    int<lower=1> p;
    vector[m] y;
    matrix[m, p] X;
    real D[m];
    int<lower=1,upper=m> i; // index for pi_i(A)
    real int_put; // input for integration
  }

  parameters {
    vector[p] beta;
    real<lower=0> A;
    vector[m] theta;;
  }

  model {
    // Prior for theta_i | beta, A
    theta ~ normal(X * beta, sqrt(A));
  
    // Likelihood: y_i | theta_i
    y ~ normal(theta, sqrt(D));
  
    // Custom prior for A using integrate_1d
    {
      vector[m] denom;
      for(j in 1:m)
        denom[j] = A + D[j];
    
      target += log(sum(inv(square(denom))));       //  ??? 1 / (A + D_u)^2
      target += 2 * log(A + D[i]);             // (A + D_i)^2
      target += log(A);                        // A
      target += -int_put;                // from the actual integration
    }
  }
  "
  
  #==============================================================================
  func_AS <- function(data_theta_y, i_val){
    # print(paste0("start of area ",i_val))
    int_func <- function(Ao){
      yy <- data_theta_y$y; Di <- D_val; D <- Di[i_val]; X <- X_mat
      num = yy[i_val] - (t(X[i_val,]) %*% 
                           ginv(t(X) %*% diag(1/(Ao+Di)) %*% X) %*% 
                           t(X)%*% diag(1/(Ao+Di)) %*% yy)
      return((num/(Ao + D))^2)
    }
    int_val <- integrate(Vectorize(int_func), 
                         lower = 0, upper = 1000, # lower = 0.005
                         rel.tol = 1e-6, abs.tol = 1e-8)$value # adjust the range of integral 
    
    options(buildtools.check=function(action) TRUE)
    # input data for area dependent prior
    data.stan <- list(m = m_real, p = p_val, i = i_val, 
                      int_put = int_val,
                      y = data_theta_y$y, X = X_mat, D = D_val)
    
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
    
    MC.summary$AS_lower <-result_AS$AS_lower[i_val]
    MC.summary$AS_upper <- result_AS$AS_upper[i_val]
    
    theta_sample <- as.vector(ex.post$theta[,i_val])
    
    MC.summary$chk_YL <- mean(ifelse(theta_sample >= MC.summary$YL_lower & 
                                       theta_sample <= MC.summary$YL_upper, 1,0))
    MC.summary$chk_AS <- mean(ifelse(theta_sample >= MC.summary$AS_lower & 
                                       theta_sample <= MC.summary$AS_upper, 1,0))
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
  AS_cov = rowMeans(matrix(unlist(lapply(1:m_real, function(i){
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_comb[[mc_itr]][[i]][[1]]$chk_AS))
  } )), byrow = T, ncol = length(M_MC))),
  YL_cov = rowMeans(matrix(unlist(lapply(1:m_real, function(i){
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_comb[[mc_itr]][[i]][[1]]$chk_YL))
  } )), byrow = T, ncol = length(M_MC))),
  med_AS = unlist(apply(matrix(unlist(lapply(1:m_real, function(i){
    unlist(lapply(M_MC, function(mc_itr) mcmc_result_comb[[mc_itr]][[i]][[1]]$chk_AS))
  } )), byrow = T, ncol = length(M_MC)),1,median)),
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

save.image("sim_unbal_comb_post_cov_S12_m15_1_100.RData")
write.csv(post_cov, "sim_unbal_comb_post_cov_S12_m30_1_40.csv", row.names = F)




