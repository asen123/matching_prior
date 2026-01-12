#==============================================================================
#  Functions to create 95% confidence interval of YL 2014
#==============================================================================
# function to return the adjusted likelihood function
YL_func <- function(m,A,Di,y,X,i){
  I <- diag(m)
  V <- diag(A+Di)
  VV1 <- diag(1/(A+Di))
  P <- VV1%*%(I-X%*%ginv(t(X)%*%VV1%*%X)%*%t(X)%*%VV1)
  D <- Di[i]
  
  int_func <- function(Ao){
    sum(1/(Ao+Di)^2)*
      (t(X[i,])%*%ginv(t(X) %*% diag(1/(Ao+Di)) %*% X) %*% X[i,])
  }
  
  int_val <- integrate(Vectorize(int_func), 
                       lower = 0, upper = 100^2, 
                       rel.tol = 1e-6, abs.tol = 1e-8)$value
  # log likelihood function to maximize
  l <- (det(V)^(-1/2))*(det(t(X)%*% VV1 %*%X)^(-1/2))*
    exp(-1/2*t(y)%*%P%*%y)
  # adjustment term
  h <- (A^((1+z^2)/4))*((A+D)^((7-z^2)/4))*exp(1/2*int_val)
  return(l*h)
}
#==============================================================================
# function to return the REML likelihood
REML_func <- function(m,A,Di,y,X){
  I <- diag(m)
  V <- diag(A+Di)
  VV1 <- diag(1/(A+Di))
  P <- VV1%*%(I-X%*%ginv(t(X)%*%VV1%*%X)%*%t(X)%*%VV1)
  # log likelihood function to maximize
  l <- (det(V)^(-1/2))*(det(t(X)%*% VV1 %*%X)^(-1/2))*
    exp(-1/2*t(y)%*%P%*%y)
  return(l)
}
#==============================================================================
# function for g1
g1 <- function(m,A,Di){
  g1 <- A*Di/(A+Di)
  return(matrix(g1))
}
#================================================================================
# function for EBLUP which uses previous function A.estimators and g1
EBLUP_func_YL <- function(m,Di,X,yy){
  # m = m_real; Di = baseball$vardir; X = X_mat; yy = data_theta_y$y
  # Find REML
  Ahat_REML <- optimize(function(A)
  { return(REML_func(m = m, A = A, Di = Di, y = yy, X = X))}, 
  interval =  c(0,100), maximum = T)$maximum
  # Find YL
  Ahat_YL <- NULL
  for(i_val in 1:m){
    Ahat_YL[i_val] <- optimize(function(A)
    { return(YL_func(m = m, A = A, Di = Di, y = yy, X = X, i=i_val))}, 
    interval =  c(0,100), maximum = T)$maximum
  }
  # combine REML, YL
  A.est <- c(Ahat_REML,Ahat_YL)
  
  # Bi values
  B <- matrix(0,m,2)
  B[,1] <- Di/(A.est[1]+Di) # B for REML
  B[,2] <- Di/(A.est[2:(m+1)]+Di) # B for AdjML
  
  x.beta.hat <- matrix(0,m,2)
  x.beta.hat[,1] <- X%*%ginv(t(X)%*%diag(B[,1]/Di)%*%X)%*%t(X)%*%diag(B[,1]/Di)%*%yy
  x.beta.hat[,2] <- X%*%ginv(t(X)%*%diag(B[,2]/Di)%*%X)%*%t(X)%*%diag(B[,2]/Di)%*%yy
  
  Eblup <- matrix(0,m,2)  
  Eblup[,1] <- yy-B[,1]*(yy-x.beta.hat[,1]) # REML eblup
  Eblup[,2] <- yy-B[,2]*(yy-x.beta.hat[,2]) # AdjML eblup
  # table of intervals (REML and YL)
  int_summ <- data.frame(A_hat_REML = rep(A.est[1],m),
                         A_hat_YL = A.est[2:(m+1)],
                         EBLUP_REML = Eblup[,1],
                         EBLUP_YL = Eblup[,2],
                         REML_lower = Eblup[,1] - z*sqrt(g1(m,A.est[1],Di)),
                         REML_upper = Eblup[,1] + z*sqrt(g1(m,A.est[1],Di)),
                         YL_lower = Eblup[,2] - z*sqrt(g1(m,A.est[2:(m+1)],Di)),
                         YL_upper = Eblup[,2] + z*sqrt(g1(m,A.est[2:(m+1)],Di)),
                         true_val = data_theta_y$theta,
                         YL_length = 2*z*sqrt(g1(m,A.est[2:(m+1)],Di)))
  # check to see if interval contains true theta
  int_summ$REML_chk <- ifelse(int_summ$true_val >= int_summ$REML_lower & 
                                int_summ$true_val <= int_summ$REML_upper,1,0)
  
  int_summ$YL_chk <- ifelse(int_summ$true_val >= int_summ$YL_lower & 
                              int_summ$true_val <= int_summ$YL_upper,1,0)
  return(int_summ)  
}
