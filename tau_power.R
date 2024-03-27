rm(list = ls())

library(parallel)
library(dplyr)
library(MASS)

numCores <- detectCores() - 1
cl <- makeCluster(numCores)

tau_power_score_fun <- function(num) {
  
  library(dplyr)
  library(MASS)
  
  
  RBNBZI <- function(n, b10, b11, b20, b21, gam10, gam11, w, tau1, tau2, grid = 20) {
    
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
  
  logL_H0 <- function(param) {
    beta1 <- param[c(1, 2)]
    beta2 <- param[c(3, 4)]
    gamma <- param[c(5, 6)]
    w <- param[7]
    
    d <- 1 - exp(-1)
    y <- sample_data$y; x1 <- sample_data$x1; x2 <- sample_data$x2; z <- sample_data$z
    
    mu1 <- c(exp(x1 %*% beta1)); mu2 <- c(exp(x2 %*% beta2))
    c1 <- exp(-d * mu1); c2 <- exp(-d * mu2)
    
    phi <- c(exp(z %*% gamma) / (1 + exp(z %*% gamma)))
    zg <- c(exp(z %*% gamma))
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    v1 <- sum(log(exp(z %*% gamma) + dpois(0, lambda = mu1) * dpois(0, lambda = mu2) * (1 + w * (1 - c1) * (1 - c2))) * ind)
    v21 <- (dpois(x = y[, 1], lambda = mu1, log = T) + dpois(x = y[, 2], lambda = mu2, log = T)) * !ind
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
  iteration = 1000
  
  result <- matrix(0, nrow = 1000, ncol = 3)
  
  n_grid <- c(200, 500)
  tau1_grid <- c(1e-10)
  tau2_grid <- c(1e-10, 0.025, 0.05, 0.1, 0.2, 0.4)
  phi_grid <- c(-0.997, 0.352)
  w_grid <- c(-1, 0, 1)
  
  param_grid <- expand.grid(w_grid, phi_grid, tau2_grid, tau1_grid, n_grid)
  
  for(i in 1:iteration) {
    set.seed(i)
    sample_data <- RBNBZI(n = param_grid[num, 5], .2, .4, .4, .8, -1.2, param_grid[num, 2], 
                          w = param_grid[num, 1], tau1 = param_grid[num, 4], tau2 = param_grid[num, 3])
    
    one <- try(optim(par = rep(0, 7), fn = logL_H0, lower = c(rep(-5, 6), -10), method = "L-BFGS-B"), silent = TRUE)
    two <- try(optim(par = rep(0, 9), fn = logL_H1, lower = c(rep(-5, 6), 1e-10, 1e-10, -10), method = "L-BFGS-B"), silent = TRUE)
    
    if (class(one) == "try-error" | class(two) == "try-error") {
      next
    } 
    
    param <- one$par
    
    ## LR-test
    L0val <- one$value
    L1val <- two$value
    LR <- 2 * (L0val - L1val)
    
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
    
    score <- colSums(gr_df)
    King_score <- colSums(gr_df)
    
    score[1:2] <- ifelse(score[1:2] <= 0, 0, score[1:2])
    score_zero <- t(score) %*% MASS::ginv(I) %*% score
    
    ## King & Wu test
    J <- MASS::ginv(I)[1:2, 1:2]
    ell <- rep(1, 2)
    JJ <- sqrt(c(t(ell) %*% MASS::ginv(J) %*% ell))
    King <- sum(King_score[1:2]) / JJ
    
    ## Result-Store
    
    result[i, 1] <- LR
    result[i, 2] <- score_zero
    result[i, 3] <- King
  }
  return(result)
}

tau_power_score <- parLapply(cl, X = 1:72, tau_power_score_fun)
stopCluster(cl)


#####################################################

load("C:/Users/jieun/Dropbox/BZINB/zero-zero/graph/tau_power_score.RData")

n_grid <- c(200, 500)
tau1_grid <- c(1e-10)
tau2_grid <- c(1e-10, 0.025, 0.05, 0.1, 0.2, 0.4)
phi_grid <- c(-0.997, 0.352)
pp_grid <- c(0.1, 0.4)
w_grid <- c(-1, 0, 1)

param_grid <- expand.grid(w_grid, phi_grid, tau2_grid, tau1_grid, n_grid)

qval = c(0.05)
upp = qchisq(qval, 2, lower.tail = F)/4 + qchisq(qval, 1, lower.tail = F)/2
kupp <- qnorm(1-qval)


# Open PNG file for writing
tiff('tau_power.tif', units = "px", res = 300, width = 2600, height = 1500)

# Set the outer margin
par(oma = c(0, 0, 0, 0))

# Set the inner margin
par(mar = c(5, 4, 4, 2) + 0.1)


par(mfrow = c(2, 3))

kk = 0
pn = 0
for(pp in phi_grid){
  pn = pn + 1
  for(ww in w_grid){
    kk = kk + 1
    
    out <- data.frame(n = rep(0, 12), tau2 = rep(0, 12), LR = rep(0, 12), score = rep(0, 12))
    
    # Score_mixed 결과
    k = 0
    for(ii in 1:72){
      if(param_grid[ii, 2] == pp & param_grid[ii, 1] == ww){  ## tau1 조정
        cat("n =", param_grid[ii, 5], "tau1 =", param_grid[ii, 4], "tau2 =", param_grid[ii, 3], "w =", param_grid[ii, 1], '\n')
        k = k + 1
        out$n[k] <- param_grid[ii, 5]
        out$tau2[k] <- param_grid[ii, 3]
        out$LR[k] <- round(mean(tau_power_score[[ii]][, 1] >= upp, na.rm = T), 3)
        out$score[k] <- round(mean(tau_power_score[[ii]][, 3] >= kupp, na.rm = T), 3)
      }
    }
    
    # if(ww == -1.5){
    #   out_full <- out  
    # } else{
    #   out_full <- rbind(out_full, out)
    # }
    print(out$score_zero[1])
    
    # 그림 그리기
    # n=200, LR
    plot(out$tau2[out$n == 500], out$LR[out$n == 500], lty = 1, type = 'l', col = 'darkgray',
         xlab = expression(tau[2]), ylab = "estimated power", main = bquote(phi == .(pp_grid[pn])  ~ ", " ~ w == .(ifelse(ww == 0, 0, ww))), ylim = c(0, 1))
    points(out$tau2[out$n == 500], out$LR[out$n == 500], lty = 1, lwd = 1, pch = 15, col = 'darkgray')
    
    # n=100, LR
    points(out$tau2[out$n == 200], out$LR[out$n == 200], lty = 2, type = 'l', col = 'darkgray',
           xlab = expression(tau[2]), ylab = "estimated power")
    points(out$tau2[out$n == 200], out$LR[out$n == 200], lty = 2, lwd = 1, pch = 15, col = 'darkgray')
    
    # n=200, score
    points(out$tau2[out$n == 500], out$score[out$n == 500], lty = 1, type = 'l', col = 'black',
           xlab = expression(tau[2]), ylab = "estimated power")
    points(out$tau2[out$n == 500], out$score[out$n == 500], lwd = 1, pch = 16, col = 'black')
    
    # n=100, score
    points(out$tau2[out$n == 200], out$score[out$n == 200], lty = 2, type = 'l', col = 'black',
           xlab = expression(tau[2]), ylab = "estimated power")
    points(out$tau2[out$n == 200], out$score[out$n == 200], lty = 2, lwd = 1, pch = 16, col = 'black')
    
    abline(h = 0.05)
    
    if(kk == 6)
      legend("bottomright", legend = c('n=500, score', 'n=200, score', 'n=500, LR', 'n=200, LR'),
             col = c("black", "black", "darkgray","darkgray"), lty = c(1, 3, 1, 3), lwd = rep(1, 4),
             pch = c(16, 16, 15, 15))
    
  }
}
dev.off()
