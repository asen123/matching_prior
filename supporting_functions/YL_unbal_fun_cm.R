#==============================================================================
#  Functions to create 95% confidence interval of YL 2014
#==============================================================================
# function to return the adjusted likelihood function
YL_func <- function(m,A,Di,y,i){
  I <- diag(m)
  V <- diag(A+Di)
  VV1 <- diag(1/(A+Di))
  D <- Di[i]
  
  # log likelihood function to maximize
  l <- (det(V)^(-1/2)) * exp(-1/2*t(y)%*% VV1 %*%y)
  # adjustment term
  h <- (A^((1+z^2)/4))*((A+D)^((7-z^2)/4))
  return(l*h)
}
#==============================================================================
# function to return the REML likelihood
REML_func <- function(m,A,Di,y){
  I <- diag(m)
  V <- diag(A+Di)
  VV1 <- diag(1/(A+Di))
  # log likelihood function to maximize
  l <- (det(V)^(-1/2))*exp(-1/2*t(y)%*% VV1 %*%y)
  return(l)
}
#==============================================================================
# function for g1
g1 <- function(m,A,Di){
  g1 <- A*Di/(A+Di)
  return(matrix(g1))
}
#================================================================================
# function for BP which uses previous function A.estimators and g1
BP_func_YL <- function(m,Di,yy){
  # m = m_real; Di = baseball$vardir; X = X_mat; yy = data_theta_y$y
  # Find REML
  Ahat_REML <- optimize(function(A)
  { return(REML_func(m = m, A = A, Di = Di, y = yy))}, 
  interval =  c(0,100), maximum = T)$maximum
  # Find YL
  Ahat_YL <- NULL
  for(i_val in 1:m){
    Ahat_YL[i_val] <- optimize(function(A)
    { return(YL_func(m = m, A = A, Di = Di, y = yy, i=i_val))}, 
    interval =  c(0,100), maximum = T)$maximum
  }
  # combine REML, YL
  A.est <- c(Ahat_REML,Ahat_YL)
  
  # Bi values
  B <- matrix(0,m,2)
  B[,1] <- Di/(A.est[1]+Di) # B for REML
  B[,2] <- Di/(A.est[2:(m+1)]+Di) # B for AdjML
  
  BP <- matrix(0,m,2)  
  BP[,1] <- (1- B[,1])*yy # ML BP
  BP[,2] <- (1- B[,2])*yy # AdjML BP
  
  # table of intervals (REML and YL)
  int_summ <- data.frame(A_hat_REML = rep(A.est[1],m),
                         A_hat_YL = A.est[2:(m+1)],
                         BP_REML = BP[,1],
                         BP_YL = BP[,2],
                         REML_lower = BP[,1] - z*sqrt(g1(m,A.est[1],Di)),
                         REML_upper = BP[,1] + z*sqrt(g1(m,A.est[1],Di)),
                         YL_lower = BP[,2] - z*sqrt(g1(m,A.est[2:(m+1)],Di)),
                         YL_upper = BP[,2] + z*sqrt(g1(m,A.est[2:(m+1)],Di)),
                         true_val = data_theta_y$theta,
                         YL_length = 2*z*sqrt(g1(m,A.est[2:(m+1)],Di)))
  # check to see if interval contains true theta
  int_summ$REML_chk <- ifelse(int_summ$true_val >= int_summ$REML_lower & 
                                int_summ$true_val <= int_summ$REML_upper,1,0)
  
  int_summ$YL_chk <- ifelse(int_summ$true_val >= int_summ$YL_lower & 
                              int_summ$true_val <= int_summ$YL_upper,1,0)
  return(int_summ)  
}
