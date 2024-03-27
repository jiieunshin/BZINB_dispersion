stopCluster(cl)
rm(list = ls()); gc(reset = T)

library(parallel)
numCores <- detectCores() - 1

cl <- makeCluster(numCores)

cond_tau_null <- function(num) {
  library(dplyr)
  library(MASS)
  
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
  
  iter <- 1000
  B <- 1000
  
  n_grid <- c(200, 500)
  tau_grid <- c(0.2, 0.4, 0.6, 0.8, 1)
  w_grid <- c(-1, 0, 1)
  
  param_grid <- expand.grid(w_grid, tau_grid, n_grid)
  
  result <- vector(length = iter)
  boot_cric <- matrix(0, nrow = iter, ncol = 4)
  
  for (i in 1:iter) {
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
    
    if(result[i] < 0.6){
      boot_cric[i, ] <- c(0, 0, 0, 0)
    } else if(result[i] > 2.5){
      boot_cric[i, ] <- c(1, 1, 1, 1)
    } else{
      ###############
      ## Bootstrap step
    K = nrow(y)
    grid = 50
    boot_b <- vector(length = B)
    
    for(b in 1:B) {
      
      boot_id = sample(1:param_grid[num, 3], param_grid[num, 3], replace = TRUE)
      yB = y[boot_id, ]
      
      x1B <- sample_data$x1[boot_id, ]
      x2B <- sample_data$x2[boot_id, ]
      zB <- sample_data$z[boot_id, ]
      indB <- (yB[, 1] == 0 & yB[, 2] == 0)
      
      logL_H0B <- function(paramB) {
        
        beta1B <- paramB[c(1, 2)]
        beta2B <- paramB[c(3, 4)]
        gammaB <- paramB[c(5, 6)]
        tau2B <- paramB[7]
        wB <- paramB[8]
        
        d <- 1 - exp(-1)
        
        mu1B <- c(exp(x1B %*% beta1B)); mu2B <- c(exp(x2B %*% beta2B))
        c1B <- exp(-d * mu1B);
        c2B <- (1 + d * mu2B * tau2B)^{-tau2B^{-1}}
        
        indB <- (yB[, 1] == 0 & yB[, 2] == 0)
        
        v1B <- log(exp(zB %*% gammaB) + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))[indB]
        v2B <- (dpois(x = yB[, 1], lambda = mu1B, log = T) + dnbinom(x = yB[, 2], mu = mu2B, size = 1/tau2B, log = T))[!indB]
        v3B <- (1 + wB * (exp(-yB[, 1]) - c1B) * (exp(-yB[, 2]) - c2B))
        
        v3B <- ifelse(v3B <= 0, log(1e-10), log(v3B))[!indB]
        
        logbpB <- sum(v1B) + sum(v2B) + sum(v3B) - sum(log(1+exp(zB %*% gammaB)))
        
        return(-logbpB)
      }
      
      oneB <- try(optim(par = rep(0, 8), lower = c(rep(-Inf, 6), 1e-10, -Inf), fn = logL_H0B, method = "L")$par, silent = TRUE)
      
      if (class(oneB) == "try-error") {
        next
      } else {
        
      }
      
      beta1B <- oneB[c(1, 2)]
      beta2B <- oneB[c(3, 4)]
      gammaB <- oneB[c(5, 6)]
      tau2B <- oneB[7]
      wB <- oneB[8]
      
      psiB <- exp(zB %*% gammaB)
      mu1B <- c(exp(x1B %*% beta1B)); mu2B <- c(exp(x2B %*% beta2B))
      
      c1B <- exp(-d * mu1B); c2B <- (1 + d * mu2B * tau2B)^{-tau2B^{-1}}
      
      a1B = sapply(1:nrow(yB), function(i){ 
        sv = 0
        if(yB[i, 1] > 0){
          for(v in 1:yB[i, 1]){
            sv = sv + (yB[i,1] - v)
          }
        }
        return (sv)
      })
      
      a2B = sapply(1:nrow(yB), function(i){ 
        sv = 0
        if(yB[i, 2] > 0){
          for(v in 1:yB[i, 2]){
            sv = sv + (yB[i,2] - v) / (1 + tau2B * yB[i, 2] - tau2B * v)
          }
        }
        return (sv)
      })
      
      dgam0B <- (psiB / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB - psiB/(1 + psiB)
      
      dgam1B <- (psiB *  zB[, 2] / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB - psiB * zB[, 2] /(1 + psiB)
      
      dbeta10B <- (((exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * mu1B * (wB * d * (1 - c2B) * c1B - (1 + wB * (1 - c1B) * (1 - c2B)))) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB) + 
        (((yB[, 1] - mu1B) + wB * d * mu1B * (exp(-yB[, 2]) - c2B) * c1B / (1 + wB * (exp(-yB[, 1]) - c1B) * (exp(-yB[, 2]) - c2B))) * !indB)
      dbeta11B <- (((exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * mu1B * x1B[, 2] * (wB * d * (1 - c2B) * c1B - (1 + wB * (1 - c1B) * (1 - c2B)))) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB) + 
        (((yB[, 1] - mu1B) * x1B[, 2] + wB * d * mu1B * x1B[, 2] * (exp(-yB[, 2]) - c2B) * c1B / (1 + wB * (exp(-yB[, 1]) - c1B) * (exp(-yB[, 2]) - c2B))) * !indB)
      dbeta20B <- ((exp(-mu1B) * mu2B * (1 + tau2B * mu2B)^{-1/tau2B} * (wB * d * (1 - c1B) * (1 + d * tau2B * mu2B)^{-1/tau2B - 1} - (1 + tau2B * mu2B)^{-1} * (1 + wB * (1 - c1B) * (1 - c2B))) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 -c2B)))) * indB) + 
        (((yB[, 2] - mu2B) / (1 + tau2B * mu2B) + wB * d * mu2B * (exp(-yB[, 1]) - c1B) * (1 + d * tau2B * mu2B)^{-1/tau2B - 1} / (1 + wB * (exp(-yB[, 1])- c1B) * (exp(-yB[, 2]) - c2B))) * !indB)
      dbeta21B <- ((exp(-mu1B) * mu2B * x2B[, 2] * (1 + tau2B * mu2B)^{-1/tau2B} * (wB * d * (1 - c1B) * (1 + d * tau2B * mu2B)^{-1/tau2B - 1} - (1 + tau2B * mu2B)^{-1} * (1 + wB * (1 - c1B) * (1 - c2B))) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 -c2B)))) * indB) + 
        (((yB[, 2] - mu2B) * x2B[, 2] / (1 + tau2B * mu2B) + wB * d * mu2B * x2B[, 2] * (exp(-yB[, 1]) - c1B) * (1 + d * tau2B * mu2B)^{-1/tau2B - 1} / (1 + wB * (exp(-yB[, 1])- c1B) * (exp(-yB[, 2]) - c2B))) * !indB)
      
      dtau1B <- (((exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)) * mu1B^2/2) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB)  -
        (((wB * c1B * (1 - c2B) * exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * d^2 * mu1B^2 / 2) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB) +
        ((a1B - yB[, 1] * mu1B + mu1B^2 / 2 - wB * c1B * (exp(-yB[, 2]) - c2B) * (d * mu1B)^2 / 2 / (1 + wB * (exp(-yB[, 1]) - c1B) * (exp(-yB[, 2]) - c2B))) * !indB)
      
      dtau2B <- ((exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)) * (tau2B^{-2} * log(1 + tau2B * mu2B) - mu2B * tau2B^{-1} / (1 + tau2B * mu2B)) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB) - 
        ((wB * (1 - c1B) * exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * c2B * (tau2B^{-2} * log(1 + d * tau2B * mu2B) - d * mu2B * tau2B^{-1} / (1 + d * tau2B * mu2B)) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB) + 
        ((a2B - yB[, 2] * mu2B / (1 + tau2B * mu2B) + tau2B^{-2} * log(1 + tau2B * mu2B) - tau2B^{-1} * mu2B / (1 + tau2B * mu2B)) * !indB) -
        ((wB * (exp(-yB[, 1]) - c1B) * c2B * (tau2B^{-2} * log(1 + d * tau2B * mu2B) - d * mu2B * tau2B^{-1} / (1 + d * tau2B * mu2B)) / (1 + wB * (exp(-yB[, 1]) - c1B) * (exp(-yB[, 2]) - c2B))) * !indB)
      dwB <- (exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 - c1B) * (1 - c2B) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB + 
        ((exp(-yB[, 1]) - c1B) * (exp(-yB[, 2]) - c2B) / (1 + wB * (exp(-yB[, 1]) - c1B) * (exp(-yB[, 2]) - c2B))) * !indB
      
      grad_dfB <- data.frame(dgam0B, dgam1B, dbeta10B, dbeta11B, dbeta20B, dbeta21B, dtau2B, dwB, dtau1B)
      
      scoreB <- colSums(grad_dfB)
      
      IB <- matrix(c(colSums(grad_dfB * dgam0B),
                     colSums(grad_dfB * dgam1B),
                     colSums(grad_dfB * dbeta10B),
                     colSums(grad_dfB * dbeta11B),
                     colSums(grad_dfB * dbeta20B),
                     colSums(grad_dfB * dbeta21B),
                     colSums(grad_dfB * dtau2B),
                     colSums(grad_dfB * dwB),
                     colSums(grad_dfB * dtau1B)), 9, 9) 
      
      
      I11B <- IB[1:8, 1:8]
      I22B <- diag(IB)[9]
      I12B <- IB[1:8, 9]
      
      boot_b[b] <- (scoreB[9] - score[9]) / sqrt(I22B - t(I12B) %*% MASS::ginv(I11B) %*% (I12B))
    }
    
    B_ctic <- quantile(boot_b, c(0.8, 0.9, 0.95, 0.99))
    boot_cric[i, ] <- ifelse(result[i] > B_ctic, 1, 0)
  } 
  
    if (i == iter + 1) {
      break
    }
  }
  return(list(result, boot_cric))
}

tau_cond_boot <- parLapply(cl, X = 1:30, cond_tau_null)
tau_cond_boot


save(result_list, file = "C:/Users/UOS/Desktop/장동민/result.Rdata")


## 결과


n_grid <- c(200, 500)
tau_grid <- c(0.2, 0.4, 0.6, 0.8, 1)
w_grid <- c(-1, 0, 1)

param_grid <- expand.grid(w_grid, tau_grid, n_grid)

qval <- c(0.99, 0.95, 0.9, 0.8)
upp = qnorm(qval)

qnorm(0.90)
grid_sets


out = list()
# 스코어 결과
for(ii in 1:30){
  cat("\n n =", param_grid[ii, 3], "gam2 =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  aaa = round(colMeans(tau_cond_boot[[ii]][[2]], na.rm = T), 3)
  cat("estimated prob =", aaa[4], aaa[3], aaa[2], aaa[1], '\n')
}

