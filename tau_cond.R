rm(list = ls())

library(parallel)
library(dplyr)
library(MASS)

# numCores <- detectCores() - 1
cl <- makeCluster(40)
result_list = list()
boot_list = list()

please <- function(num) {
  # for(num in 2:15){
  library(dplyr)
  library(MASS)
  
  # n = param_grid[num, 3]; b10 = .2; b11 = .4; b20 = .4; b21 = .8; gam10 = -1.2; gam11 = 2.4; w = param_grid[num, 2]; tau1 = 0.8; tau2 = param_grid[num, 1]
  
  RBNBZI <- function(n, b10, b11, b20, b21, gam10, gam11, w, tau1, tau2, grid = 50) {
    
    find_idx <- function(mat) {
      if(sum(mat) == 0){
        return(c(grid-1, grid-1))
      } else{
        for (row in 1:grid) {
          for (col in 1:grid) {
            if (mat[row, col] == TRUE) {
              return(c(row - 1, col - 1))
            }
          }
        }
      }
    }
    
    d <- 1 - exp(-1)
    
    x1 <- cbind(1, runif(n))
    x2 <- cbind(1, runif(n))
    z <- cbind(1, runif(n))
    
    mu1 <- exp(x1 %*% c(b10, b11));
    mu2 <- exp(x2 %*% c(b20, b21));
    phi <- exp(z %*% c(gam10, gam11)) / (1 + exp(z %*% c(gam10, gam11)))
    
    c1 <- (1 + d * mu1 * tau1)^{-1/tau1}
    c2 <- (1 + d * mu2 * tau2)^{-1/tau2}
    
    y <- matrix(0, nrow = n, ncol = 2)
    
    for (i in 1:n) {
      PX <- cbind(dnbinom(0:(grid - 1), mu = mu1[i], size = 1/tau1), dnbinom(0:(grid - 1), mu = mu2[i], size = 1/tau2))
      PTMPY <- cbind(exp(-(0:(grid - 1))) - (1 + d * mu1[i] * tau1)^{-1/tau1}, exp(-(0:(grid - 1))) - (1 + d * mu2[i] * tau2)^{-1/tau2})
      
      TXY = matrix(PX[, 1], nrow = grid, ncol = 1) %*% matrix(PX[, 2], nrow = 1, ncol = grid)  * 
        (1 + w * matrix(PTMPY[, 1], nrow = grid, ncol = 1) %*% matrix(PTMPY[, 2], nrow = 1, ncol = grid))
      
      SXY = matrix(cumsum(TXY %>% t() %>% as.vector()), nrow = grid, ncol = grid, byrow = TRUE)
      
      RAN_N <- runif(1)
      
      if (RAN_N < phi[i]) {
        y[i, ] <- c(0, 0)
      } else {
        RAN_A <- runif(1)
        mat = RAN_A < SXY
        find_idx(mat)
        y[i, ] <- find_idx(mat)
      }
    }
    return(list(y = y, x1 = x1, x2 = x2, z = z, c1 = c1, c2 = c2))
  }
  
  logL_H0 <- function(param) {
    beta1 <- param[c(1, 2)]
    beta2 <- param[c(3, 4)]
    gamma <- param[c(5, 6)]
    tau <- c(1e-10, param[7])
    w <- param[8]
    
    d <- 1 - exp(-1)
    y <- sample_data$y; x1 <- sample_data$x1; x2 <- sample_data$x2; z <- sample_data$z
    
    mu1 <- c(exp(x1 %*% beta1)); mu2 <- c(exp(x2 %*% beta2))
    c1 <- exp(-d * mu1); c2 <- (1 + d * mu2 * tau[2])^{-tau[2]^{-1}}
    
    phi <- c(exp(z %*% gamma) / (1 + exp(z %*% gamma)))
    zg <- c(exp(z %*% gamma))
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    v1 <- sum(log(exp(z %*% gamma) + dnbinom(0, mu = mu2, size = 1/tau[2]) * dpois(0, mu1) * (1 + w * (1 - c1) * (1 - c2))) * ind)
    
    v21 <- ( dpois(y[, 1], mu1, log = T) + dnbinom(x = y[, 2], mu = mu2, size = 1/tau[2], log = T)) * !ind
    v22 <- c((1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2)) * !ind)
    v22 <- ifelse(v22 <= 0, log(1e-15), log(v22)) * !ind
    v2 <- (v21 + v22)*!ind
    logbp <- sum(v1) + sum(v2) - sum(log(1+exp(z %*% gamma)))
    
    return(-logbp)
  }
  
  logL_H1 <- function(param) {
    beta1 <- param[c(1, 2)]
    beta2 <- param[c(3, 4)]
    gamma <- param[c(5, 6)]
    tau <- param[c(7, 8)]
    w <- param[9]
    
    d <- 1 - exp(-1)
    y <- sample_data$y; x1 <- sample_data$x1; x2 <- sample_data$x2; z <- sample_data$z
    
    mu1 <- c(exp(x1 %*% beta1)); mu2 <- c(exp(x2 %*% beta2))
    c1 <- (1 + d * mu1 * tau[1])^{-1/tau[1]}; c2 <- (1 + d * mu2 * tau[2])^{-1/tau[2]}
    
    phi <- c(exp(z %*% gamma) / (1 + exp(z %*% gamma)))
    zg <- c(exp(z %*% gamma))
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    v1 <- sum(log(exp(z %*% gamma) + dnbinom(0, mu = mu1, size = 1/tau[1]) * dnbinom(0, mu = mu2, size = 1/tau[2]) * (1 + w * (1 - c1) * (1 - c2))) * ind)
    
    v21 <- (dnbinom(x = y[, 1], mu = mu1, size = 1/tau[1], log = T) + dnbinom(x = y[, 2], mu = mu2, size = 1/tau[2], log = T)) * !ind
    v22 <- c((1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2)) * !ind)
    v22 <- ifelse(v22 <= 0, log(1e-15), log(v22)) * !ind
    v2 <- (v21 + v22)*!ind
    logbp <- sum(v1) + sum(v2) - sum(log(1+exp(z %*% gamma)))
    
    return(-logbp)
  }
  
  ## Iteration
  # i = 1
  iteration = 100000
  result <- c()
  
  n_grid <- c(200, 500)
  tau_grid <- c(0.2, 0.4, 0.6, 0.8, 1)
  w_grid <- c(-1, 0, 1, 2)
  
  param_grid <- expand.grid(w_grid, tau_grid, n_grid)
  param_est <- matrix(0, iteration, 8)
  for(i in 1:iteration) {
    set.seed(i)
    sample_data <- RBNBZI(n = param_grid[num, 3], b10 = .2, b11 = .4, b20 = .4, b21 = .8, gam10 = -1.2, gam11 = 2.4, w = param_grid[num, 1], tau1 = 1e-10, tau2 = param_grid[num, 2] )
    
    one <- try(optim(par = rep(0, 8), lower = c(rep(-Inf, 6), 1e-10, -Inf), fn = logL_H0, method = "L")$par, silent = TRUE)
    
    if(class(one) == "try-error") {
      next
    }
    
    param <- one
    
    d <- 1 - exp(-1)
    
    y <- sample_data$y; x1 <- sample_data$x1; x2 <- sample_data$x2; z <- sample_data$z
    
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    beta1 <- param[c(1, 2)]
    beta2 <- param[c(3, 4)]
    gamma <- param[c(5, 6)]
    tau2 <- param[7]
    w <- param[8]
    
    psi <- exp(z %*% gamma)
    mu1 <- c(exp(x1 %*% beta1)); mu2 <- c(exp(x2 %*% beta2))
    
    c1 <- exp(-d * mu1); c2 <- (1 + d * mu2 * tau2)^{-tau2^{-1}}
    
    a1 = sapply(1:nrow(y), function(i){ 
      sv = 0
      if(y[i, 1] > 0){
        for(v in 1:y[i, 1]){
          sv = sv + (y[i,1] - v)
        }
      }
      return (sv)
    })
    
    a2 = sapply(1:nrow(y), function(i){ 
      sv = 0
      if(y[i, 2] > 0){
        for(v in 1:y[i, 2]){
          sv = sv + (y[i,2] - v) / (1 + tau2 * y[i, 2] - tau2 * v)
        }
      }
      return (sv)
    })
    
    dgam0 <- (psi / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind - psi/(1 + psi)
    
    dgam1 <- (psi *  z[, 2] / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind - psi * z[, 2] /(1 + psi)
    
    dbeta10 <- (((exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * mu1 * (w * d * (1 - c2) * c1 - (1 + w * (1 - c1) * (1 - c2)))) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind) + 
      (((y[, 1] - mu1) + w * d * mu1 * (exp(-y[, 2]) - c2) * c1 / (1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2))) * !ind)
    dbeta11 <- (((exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * mu1 * x1[, 2] * (w * d * (1 - c2) * c1 - (1 + w * (1 - c1) * (1 - c2)))) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind) + 
      (((y[, 1] - mu1) * x1[, 2] + w * d * mu1 * x1[, 2] * (exp(-y[, 2]) - c2) * c1 / (1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2))) * !ind)
    dbeta20 <- ((exp(-mu1) * mu2 * (1 + tau2 * mu2)^{-1/tau2} * (w * d * (1 - c1) * (1 + d * tau2 * mu2)^{-1/tau2 - 1} - (1 + tau2 * mu2)^{-1} * (1 + w * (1 - c1) * (1 - c2))) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 -c2)))) * ind) + 
      (((y[, 2] - mu2) / (1 + tau2 * mu2) + w * d * mu2 * (exp(-y[, 1]) - c1) * (1 + d * tau2 * mu2)^{-1/tau2 - 1} / (1 + w * (exp(-y[, 1])- c1) * (exp(-y[, 2]) - c2))) * !ind)
    dbeta21 <- ((exp(-mu1) * mu2 * x2[, 2] * (1 + tau2 * mu2)^{-1/tau2} * (w * d * (1 - c1) * (1 + d * tau2 * mu2)^{-1/tau2 - 1} - (1 + tau2 * mu2)^{-1} * (1 + w * (1 - c1) * (1 - c2))) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 -c2)))) * ind) + 
      (((y[, 2] - mu2) * x2[, 2] / (1 + tau2 * mu2) + w * d * mu2 * x2[, 2] * (exp(-y[, 1]) - c1) * (1 + d * tau2 * mu2)^{-1/tau2 - 1} / (1 + w * (exp(-y[, 1])- c1) * (exp(-y[, 2]) - c2))) * !ind)
    
    dtau1 <- (((exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)) * mu1^2/2) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind)  -
      (((w * c1 * (1 - c2) * exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * d^2 * mu1^2 / 2) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind) +
      ((a1 - y[, 1] * mu1 + mu1^2 / 2 - w * c1 * (exp(-y[, 2]) - c2) * (d * mu1)^2 / 2 / (1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2))) * !ind)
    
    dtau2 <- ((exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)) * (tau2^{-2} * log(1 + tau2 * mu2) - mu2 * tau2^{-1} / (1 + tau2 * mu2)) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind) - 
      ((w * (1 - c1) * exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * c2 * (tau2^{-2} * log(1 + d * tau2 * mu2) - d * mu2 * tau2^{-1} / (1 + d * tau2 * mu2)) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind) + 
      ((a2 - y[, 2] * mu2 / (1 + tau2 * mu2) + tau2^{-2} * log(1 + tau2 * mu2) - tau2^{-1} * mu2 / (1 + tau2 * mu2)) * !ind) -
      ((w * (exp(-y[, 1]) - c1) * c2 * (tau2^{-2} * log(1 + d * tau2 * mu2) - d * mu2 * tau2^{-1} / (1 + d * tau2 * mu2)) / (1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2))) * !ind)
    dw <- (exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 - c1) * (1 - c2) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind + 
      ((exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2) / (1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2))) * !ind
    
    grad_df <- data.frame(dgam0, dgam1, dbeta10, dbeta11, dbeta20, dbeta21, dtau2, dw, dtau1)
    score <- colSums(grad_df)
    
    I <- matrix(c(colSums(grad_df * dgam0),
                  colSums(grad_df * dgam1),
                  colSums(grad_df * dbeta10),
                  colSums(grad_df * dbeta11),
                  colSums(grad_df * dbeta20),
                  colSums(grad_df * dbeta21),
                  colSums(grad_df * dtau2),
                  colSums(grad_df * dw),
                  colSums(grad_df * dtau1)), 9, 9) 
    
    
    I11 <- I[1:8, 1:8]
    I22 <- diag(I)[9]
    I12 <- I[1:8, 9]
    
    result[i] <- score[9] / sqrt(I22 - t(I12) %*% MASS::ginv(I11) %*% (I12))
    
  }
  
  return(result)
  # return (result)
}


tau_result_list <- parLapply(cl, X = 1:40, please)
tau_result_list
stopCluster(cl)


## 결과
n_grid <- c(200, 500)
tau_grid <- c(0.2, 0.4, 0.6, 0.8, 1)
w_grid <- c(-1, 0, 1, 2)

param_grid <- expand.grid(w_grid, tau_grid, n_grid)

qval = c(0.2, 0.1, 0.05, 0.01)
upp = qchisq(qval, 2, lower.tail = F)/4 + qchisq(qval, 1, lower.tail = F)/2
upp = qchisq(qval, 2, lower.tail = F)
result_list[[1]]$result
kupp = qnorm(1-qval)

# 스코어 결과
for(ii in 1:40){
  cat("\n n =", param_grid[ii, 3], "gam2 =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  a1 = round(mean(tau_result_list[[ii]]>= kupp[4], na.rm = T), 3)
  a2 = round(mean(tau_result_list[[ii]]>= kupp[3], na.rm = T), 3)
  a3 = round(mean(tau_result_list[[ii]]>= kupp[2], na.rm = T), 3)
  a4 = round(mean(tau_result_list[[ii]]>= kupp[1], na.rm = T), 3)
  cat(a1, a2, a3, a4)
}


#histogram
id = c(16, 19, 22, 25, 28, 31, 34, 37, 40, 43)
par(mfrow=c(2,5))
for(ii in id){
  y = tau_result_list1000[[ii]][[3]]
  qqplot(qnorm(ppoints(500)), y, main = paste0("n=", param_grid[ii,3], ", tau2=", param_grid[ii,2], ", w=", param_grid[ii,1]), ylab = 'Simulated score tests', xlab = 'N(0,1) quantiles')
  qqline(y, distribution = function(p) qnorm(p), probs = c(0.001, 0.999), col = 2)
}


###################

tiff('C:/Users/jieun/Dropbox/BZINB/zero-zero/graph/tau_shift.tif', units = "px", res = 300, width = 1078*3, height = 532*3)
phi_grid <- c(0.2, 0.5)
w_grid <- c(-1, 0, 1, 2)
param_grid <- expand.grid(w_grid, phi_grid)
par(mfrow = c(2,2))
for(ii in c(1,3,5,7)){
  hist(tau_result_list1000[[ii]][[3]], breaks = 100, main = bquote(phi == .(param_grid[ii, 2]) ~ omega == .(param_grid[ii, 1])), 
       freq = F, ylim = c(0, 0.4), xlim = c(-4, 4), xlab = 'score')
  x  <- seq(from = -6, to = 6, by = 0.01)
  fx <- dnorm(x, mean = 0, sd = 1)
  lines(x, fx, type = "l", lwd = 2, col = "red")
}
dev.off()

###################

save(tau_result_list, file = "~/approx_cond_tau.RData")
load("~/tau_modify_score_n50.RData")
