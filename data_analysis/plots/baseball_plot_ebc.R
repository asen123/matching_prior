library(psych); library(MASS); library(mvtnorm); library(rstan); 
library(ggplot2); library(dplyr); library(foreach);
library(doParallel); library(tidyverse); library(parallel)
#==============================================================================
# Source YL codes
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
# baseball = read.table("/Users/aditisen/Downloads/baseball.txt")
baseball = read.table("baseball.txt")
colnames(baseball) = c("player", "est_dir", "x1", "x2", "actual")
# arc sine transformation
baseball$y <- sqrt(45)*asin(2*baseball$est_dir - 1)
baseball$theta <- sqrt(45)*asin(2*baseball$actual - 1)

# common mean
mu <- -3.275
A_real <- 0.26582 # A = prior variance
z <- 1.96 
m_real <- nrow(baseball)
data_theta_y <- list(theta = baseball$theta - mu, 
                     y = baseball$y - mu)
D_val <- rep(1,m_real) 
# Run function to generate I_YL CI
result_YL_S2 <- BP_func_YL(m = m_real, Di = D_val, yy = data_theta_y$y)
remove(data_theta_y)

# Covariate
source("AS_bal_fun.R")
source("YL_bal_fun.R")
data_theta_y <- list(theta = baseball$theta, y = baseball$y)
X_mat <- cbind(rep(1,m_real), baseball$x1, baseball$x2)
beta_ols <- solve(t(X_mat)%*%X_mat)%*%t(X_mat)%*%(data_theta_y$y)
D_val <- rep(1,m_real) 
p_val <- ncol(X_mat)

result_AS_S1 <- EBLUP_func_AS(m = m_real, Di = D_val,
                           X = X_mat, yy = data_theta_y$y)
result_YL_S1 <- EBLUP_func_YL(m = m_real, Di = D_val,
                           X = X_mat, yy = data_theta_y$y)
remove(data_theta_y)
#==============================================================================
# EB CI
#==============================================================================

EBC_baseball_S2 <- data.frame(
  area_name = baseball$player,
  true = baseball$actual,
  # YL
  YL_lower = (sin((result_YL_S2$YL_lower+mu)/sqrt(45))+1)/2,
  YL_upper = (sin((result_YL_S2$YL_upper+mu)/sqrt(45))+1)/2,
  YL_length = result_YL_S2$YL_length,
  YL_chk = result_YL_S2$YL_chk
)

EBC_baseball_S1 <- data.frame(
  area_name = baseball$player,
  true = baseball$actual,
  # AS
  AS_lower = (sin((result_AS_S1$AS_lower)/sqrt(45))+1)/2,
  AS_upper = (sin((result_AS_S1$AS_upper)/sqrt(45))+1)/2,
  AS_length = result_AS_S1$AS_length,
  AS_chk = result_AS_S1$AS_chk,
  # YL
  YL_lower = (sin((result_YL_S1$YL_lower)/sqrt(45))+1)/2,
  YL_upper = (sin((result_YL_S1$YL_upper)/sqrt(45))+1)/2,
  YL_length = result_YL_S1$YL_length,
  YL_chk = result_YL_S1$YL_chk
)

#==============================================================================

m <- m_real
sep_val <- 0.5

pdf("baseball_ebc.pdf")
plot(2*c(1:m), EBC_baseball_S1$true, type = "n", 
     ylim = c(0.1, 0.5), xlim = c(1.5, 2*m),
     xaxt = 'n', yaxt = "n",
     xlab = "", ylab = "", cex.axis = 2, cex.lab = 2)
# Draw segments first
segments(x0 = 2*c(1:m)-sep_val, y0 = EBC_baseball_S2$YL_lower, 
         x1 = 2*c(1:m)-sep_val, y1 = EBC_baseball_S2$YL_upper, 
         col = 'darkgreen', lwd = 2)
segments(x0 = 2*c(1:m), y0 = EBC_baseball_S1$AS_lower, 
         x1 = 2*c(1:m), y1 = EBC_baseball_S1$AS_upper, 
         col = 'blue', lwd = 2)
segments(x0 = 2*c(1:m)+sep_val, y0 = EBC_baseball_S1$YL_lower, 
         x1 = 2*c(1:m)+sep_val, y1 = EBC_baseball_S1$YL_upper, 
         col = 'red', lwd = 2)
# Then draw black points on top
points(2*c(1:m)-sep_val, EBC_baseball_S1$true, pch = 20, col = "black")
points(2*c(1:m), EBC_baseball_S1$true, pch = 20, col = "black")
points(2*c(1:m)+sep_val, EBC_baseball_S1$true, pch = 20, col = "black")
# Axes
axis(1, at = 2*c(1:m), labels = EBC_baseball_S1$area_name, cex.axis = 1.1, las = 2)
axis(2, at = seq(0.1, 0.5, by = 0.1), labels = TRUE, cex.axis = 2)
# Legend
legend(
  x = "topright",
  pch = 19,
  col = c("darkgreen", "blue", "red", "black"),
  legend = c(
    expression(I[i]^YL * " (M3)"),
    expression(I[i]^N * " (M4)"),
    expression(I[i]^YL * " (M4)"),
    "True"
  ),
  cex = 1.3
)
dev.off()


#==============================================================================
pdf("baseball_Ahat.pdf")
plot(1:m_real, result_AS_S1$A_hat_AS, 
     type = "n",
     xaxt = 'n', yaxt = "n",
     xlab = "", ylab = "", cex.axis = 2, cex.lab = 2,
     xlim = c(1,m), ylim = c(0.5, 0.8))
points(1:m_real, result_YL_S2$A_hat_YL, pch = 19, col = "darkgreen")
points(1:m_real, result_AS_S1$A_hat_AS,  pch = 19, col = "blue")
points(1:m_real, result_YL_S1$A_hat_YL, pch = 19, col = "red")
axis(1, at = 1:m_real, labels = baseball$player, cex.axis = 1.1, las = 2)
axis(2, at = seq(0.5,2, by = 0.1), labels = TRUE, cex.axis = 2)
legend(
  x = "bottomright",
  pch = 19,
  col = c("darkgreen","blue", "red"),
  legend = c(
    expression(tilde(A)[paste(i, ",sp")] ~ "(M3)"),
    expression(tilde(A)[i] ~ "(M4)"),
    expression(hat(A)[i] ~ "(M4)")),
  cex = 1.3
)
dev.off()
