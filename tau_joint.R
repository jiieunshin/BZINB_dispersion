rm(list = ls())

library(parallel)
library(dplyr)
library(MASS)

# numCores <- detectCores() - 1
cl <- makeCluster(75)
result_list = list()
boot_list = list()

please <- function(num) {
  # for(num in 2:15){
  library(dplyr)
  library(MASS)
  RBZI <- function(n, b10, b11, b20, b21, gam10, gam11, w, grid = 20) {
    
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
    
    c1 <- exp(-d * mu1);
    c2 <- exp(-d * mu2);
    
    y <- matrix(0, nrow = n, ncol = 2)
    
    for (i in 1:n) {
      PX <- cbind(dpois(0:(grid - 1), lambda = mu1[i]), dpois(0:(grid - 1), lambda = mu2[i]))
      PTMPY <- cbind(exp(-(0:(grid - 1))) - exp(-d * mu1[i]), exp(-(0:(grid - 1))) - exp(-d * mu2[i]))
      
      TXY = matrix(PX[, 1], nrow = grid, ncol = 1) %*% matrix(PX[, 2], nrow = 1, ncol = grid)  * 
        (1 + w * matrix(PTMPY[, 1], nrow = grid, ncol = 1) %*% matrix(PTMPY[, 2], nrow = 1, ncol = grid))
      
      SXY = matrix(cumsum(TXY %>% t() %>% as.vector()), nrow = 20, ncol = 20, byrow = TRUE)
      
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
  
  logL <- function(param) {
    beta1 <- param[c(1, 2)]
    beta2 <- param[c(3, 4)]
    gamma <- param[c(5, 6)]
    w <- param[7]
    
    d <- 1 - exp(-1)
    y <- sample_data$y; x1 <- sample_data$x1; x2 <- sample_data$x2; z <- sample_data$z
    
    mu1 <- c(exp(x1 %*% beta1)); mu2 <- c(exp(x2 %*% beta2))
    c1 <- exp(mu1 * exp(-1) - mu1); c2 <- exp(mu2 * exp(-1) - mu2)
    
    phi <- c(exp(z %*% gamma) / (1 + exp(z %*% gamma)))
    zg <- c(exp(z %*% gamma))
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    v1 <- sum(log(exp(z %*% gamma) + dpois(x = 0, lambda = mu1) * dpois(x = 0, lambda = mu2) * (1 + w * (1 - c1) * (1 - c2))) * ind)
    
    v21 <- (dpois(x = y[, 1], lambda = mu1, log = T) + dpois(x = y[, 2], lambda = mu2, log = T)) * !ind
    v22 <- c((1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2)) * !ind)
    v22 <- ifelse(v22 <= 0, log(1e-15), log(v22)) * !ind
    v2 <- (v21 + v22)*!ind
    logbp <- sum(v1) + sum(v2) - sum(log(1+exp(z %*% gamma)))
    
    return(-logbp)
  }
  
  
  ## Iteration
  # i = 1
  iteration = 4000
  result <- result2 <- modify_result <- c()
  boot_cric <- matrix(0, iteration, 4)
  
  n_grid <- c(100, 200, 500)
  g2_grid <- c(-1.994, -0.373, 0.705, 1.589, 2.4)
  w_grid <- c(-1.5, -1, 0, 1, 1.5, 2)
  
  param_grid <- expand.grid(w_grid, g2_grid, n_grid)
  param_est <- matrix(0, iteration, 7)
  for(i in 1:iteration) {
    set.seed(i)
    sample_data <- RBZI(n = param_grid[num, 3], .2, .4, .4, .6, -1.2, param_grid[num, 2], param_grid[num, 1])
    
    one <- try(optim(par = rep(0, 7), logL, method = "L-BFGS-B"), silent = TRUE)
    
    if (class(one) == "try-error") {
      next
    } else {
      param <- one$par
      param_est[i,] <- param
    }
    
    beta_i <- param[c(1, 3)]
    beta_c <- param[c(2, 4)]
    gam <- param[(c(5, 6))]
    w <- param[7]
    
    d <- 1 - exp(-1)
    
    ## Data
    y <- sample_data$y;
    x <- cbind(sample_data$x1[, 2], sample_data$x2[, 2]);
    z <- sample_data$z[,2]
    
    ## Parameter
    mu <- cbind(exp(sample_data$x1 %*% c(beta_i[1], beta_c[1])), exp(sample_data$x2 %*% c(beta_i[2], beta_c[2])))
    phi <- exp(sample_data$z %*% gam)/(1+exp(sample_data$z %*% gam))
    psi <- phi / (1 - phi)
    c <- cbind(exp(-d * mu[, 1]), exp(-d * mu[, 2]))
    D <- (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    
    ## Indicator 
    Ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    ## Gradient
    
    ## GAMMA
    dgam0 <- (exp(cbind(1, z) %*% gam) / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind -
      (exp(cbind(1, z) %*% gam) / (1 + exp(cbind(1, z) %*% gam)))
    
    dgam1 <- (z * exp(cbind(1, z) %*% gam) / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind -
      (z * exp(cbind(1, z) %*% gam) / (1 + exp(cbind(1, z) %*% gam)))
    
    
    ## BETA1
    dbeta10 <- (exp(-mu[, 1] - mu[, 2]) * mu[, 1] * (w * d * (1 - exp(-d * mu[, 2])) * exp(-d * mu[, 1]) - D) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((y[, 1] - mu[, 1]) + w * d * mu[, 1] * (exp(-y[, 2]) - exp(-d * mu[, 2])) * exp(-d * mu[, 1]) / D) * !Ind
    
    dbeta11 <- (x[, 1] * exp(-mu[, 1] - mu[, 2]) * mu[, 1] * (w * d * (1 - exp(-d * mu[, 2])) * exp(-d * mu[, 1]) - D) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((y[, 1] - mu[, 1]) * x[, 1] + w * d * mu[, 1] * x[, 1] * (exp(-y[, 2]) - exp(-d * mu[, 2])) * exp(-d * mu[, 1]) / D) * !Ind
    
    
    ## BETAA2
    dbeta20 <- (exp(-mu[, 1] - mu[, 2]) * mu[, 2] * (w * d * (1 - exp(-d * mu[, 1])) * exp(-d * mu[, 2]) - D) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((y[, 2] - mu[, 2]) + w * d * mu[, 2] * (exp(-y[, 1]) - exp(-d * mu[, 1])) * exp(-d * mu[, 2]) / D) * !Ind
    
    dbeta21 <- (x[, 2] * exp(-mu[, 1] - mu[, 2]) * mu[, 2] * (w * d * (1 - exp(-d * mu[, 1])) * exp(-d * mu[, 2]) - D) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((y[, 2] - mu[, 2]) * x[, 2] + w * d * mu[, 2] * x[, 2] * (exp(-y[, 1]) - exp(-d * mu[, 1])) * exp(-d * mu[, 2]) / D) * !Ind
    
    
    ## OMEGA
    dw = (exp(-mu[, 1] - mu[, 2]) * (1 - exp(-d * mu[, 1])) * (1 - exp(-d * mu[, 2])) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((exp(-y[, 1]) - exp(-d * mu[, 1])) * (exp(-y[, 2]) - exp(-d * mu[, 2])) / D) * !Ind
    
    
    ## TAU
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
          sv = sv + (y[i,2] - v)
        }
      }
      return (sv)
    })
    
    dtau1 <- (exp(-mu[, 1] -mu[ ,2]) * D * (mu[, 1]^2 / 2) / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind -
      (w * exp(-mu[, 1] - mu[, 2]) * exp(-d * mu[, 1]) * (1 - exp(-d * mu[, 2])) * (mu[, 1]^2 * d / 2)  / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      (a1 - mu[, 1] * y[, 1] + mu[, 1]^2/ 2 - w * (exp(-y[, 2]) - c[, 2]) * exp(-d * mu[, 1]) * (d^2 * mu[, 1]^2 / 2) / (1 + D)) * !Ind
    
    dtau2 <- (exp(-mu[, 1] -mu[ ,2]) * D * (mu[, 2]^2 / 2) / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind -
      (w * exp(-mu[, 1] - mu[, 2]) * exp(-d * mu[, 2]) * (1 - exp(-d * mu[, 1])) * (mu[, 2]^2 * d / 2)  / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      (a2 - mu[, 2] * y[, 2] + mu[, 2]^2/ 2 - w * (exp(-y[, 1]) - c[, 1]) * exp(-d * mu[, 2]) * (d^2 * mu[, 2]^2 / 2) / (1 + D)) * !Ind
    
    
    gr_df <- data.frame(dtau1, dtau2, dgam0, dgam1, dbeta10, dbeta11, dbeta20, dbeta21, dw)
    
    I <- matrix(c(colSums(gr_df * dtau1),
                  colSums(gr_df * dtau2),
                  colSums(gr_df * dgam0),
                  colSums(gr_df * dgam1),
                  colSums(gr_df * dbeta10),
                  colSums(gr_df * dbeta11),
                  colSums(gr_df * dbeta20),
                  colSums(gr_df * dbeta21),
                  colSums(gr_df * dw)), 9, 9) 
    
    EI <- I / param_grid[num, 3]
    
    
    score <- colSums(gr_df)
    
    sing <- try(ginv(I), silent = TRUE)
    sing2 <- try(ginv(EI), silent = TRUE)
    
    if (class(sing)[[1]] == "try-error") {
      next
    } else {
      
    }
    
    if (class(sing2)[[1]] == "try-error") {
      next
    } else {
      
    }
    
    result[i] <- t(score) %*% ginv(I) %*% score
    
    score2 <- score
    score2[1:2] <- ifelse(score2[1:2] < 0, 0, score[1:2])
    result2[i] <- t(score2) %*% ginv(I) %*% score2
    
    sum_score <- sum(score[1:2])
    # V <- I[1:2,1:2] / param_grid[num, 3]
    
    inv_EI <- ginv(I)[1:2,1:2]  ## J''
    ell <- rep(1, 2)
    variance <- sqrt(c(t(ell) %*% ginv(inv_EI) %*% ell))
    modify_result[i] <- sum_score/variance
    if (i %% 100 == 0) print(i)
  }

  return(list(result, result2, modify_result))
  # return (result)
}


tau_result_list1000 <- parLapply(cl, X = 1:75, please)
tau_result_list
stopCluster(cl)


## 결과
n_grid <- c(100, 200, 500)
phi_grid <- c(0.1, 0.2, 0.3, 0.4, 0.5)
w_grid <- c(-1.5, -1, 0, 1, 1.5, 2)

param_grid <- expand.grid(w_grid, phi_grid, n_grid)

qval = c(0.1, 0.05)
upp = qchisq(qval, 2, lower.tail = F)/4 + qchisq(qval, 1, lower.tail = F)/2
upp = qchisq(qval, 2, lower.tail = F)
result_list[[1]]$result
kupp = qnorm(1-qval)

# 스코어 결과
for(ii in 1:75){
  cat("n =", param_grid[ii, 3], "gam2 =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  for(qq in 1:2){
    cat("true q = ", qval[qq], "estimated prob =", round(mean(tau_result_list1000[[ii]][[3]]>= kupp[qq], na.rm = T), 3), '\n')
  }
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

save(result_list, file = "~/tau_modify_score10000.RData")
load("~/tau_modify_score_n50.RData")
