rm(list = ls())

library(parallel)
library(dplyr)
library(MASS)

numCores <- detectCores() - 1
cl <- makeCluster(numCores)
result_list = list()
boot_list = list()

please <- function(num) {
  library(dplyr)
  library(MASS)
  RBZI <- function(n, b10, b11, b20, b21, gam10, gam11, w, grid = 50) {
    
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
    
    v1 <- exp(z %*% gamma) + dpois(x = 0, lambda = mu1) * dpois(x = 0, lambda = mu2) * (1 + w * (1 - c1) * (1 - c2))
    v1 <- ifelse(v1 <= 0, 1e-15, v1)
    v1 <- sum(log(v1) * ind)
    
    v21 <- (dpois(x = y[, 1], lambda = mu1, log = T) + dpois(x = y[, 2], lambda = mu2, log = T)) * !ind
    v22 <- c((1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2)) * !ind)
    v22 <- ifelse(v22 <= 0, log(1e-15), log(v22)) * !ind
    v2 <- (v21 + v22)*!ind
    logbp <- sum(v1) + sum(v2) - sum(log(1+exp(z %*% gamma)))
    
    return(-logbp)
  }
  
  tau_gradient = function(y, x, z, mu, gam, w){
    phi <- exp(cbind(1, z) %*% gam)/(1+exp(cbind(1, z) %*% gam))
    psi <- phi / (1 - phi)
    c <- cbind(exp(-d * mu[, 1]), exp(-d * mu[, 2]))
    D <- (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    
    ## Indicator 
    Ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    # gradient
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
      (w * exp(-mu[, 1] - mu[, 2]) * exp(-d * mu[, 1]) * (1 - exp(-d * mu[, 2])) * (mu[, 1]^2 * d^2 / 2)  / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      (a1 - mu[, 1] * y[, 1] + mu[, 1]^2/ 2 - w * (exp(-y[, 2]) - c[, 2]) * exp(-d * mu[, 1]) * (d^2 * mu[, 1]^2 / 2) / (1 + D)) * !Ind
    
    dtau2 <- (exp(-mu[, 1] -mu[ ,2]) * D * (mu[, 2]^2 / 2) / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind -
      (w * exp(-mu[, 1] - mu[, 2]) * exp(-d * mu[, 2]) * (1 - exp(-d * mu[, 1])) * (mu[, 2]^2 * d^2 / 2)  / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
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
    
    return(list(G = gr_df, I = I))
  }
  
  
  ## Iteration #################################################################
  # i = 1
  iteration = 1000
  B = 1000
  result <- c()
  boot_cric <- matrix(NA, B, 4)
  
  # n_grid <- c(200, 500)
  
  n_grid <- c(200, 500)
  g2_grid <- c(-1.994, -0.373, 0.705, 1.589, 2.4)
  w_grid <- c(-1, 0, 1, 2)
  
  param_grid <- expand.grid(w_grid, g2_grid, n_grid)
  out = list()
  # for(num in 1:40){
  
  for(i in 1:iteration) {
    set.seed(i)
    sample_data <- RBZI(n = param_grid[num, 3], .2, .4, .4, .8, -1.2, param_grid[num, 2], param_grid[num, 1])
    
    one <- try(optim(par = rep(0, 7), logL, method = "L-BFGS-B"), silent = TRUE)
    
    if (class(one) == "try-error") {
      next
    } else {
      param <- one$par
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
    
    ## Gradient
    gtau_out = tau_gradient(y, x, z, mu, gam, w)
    
    if (class(ginv(gtau_out$I))[[1]] == "try-error") {
      next
    } else {
      
    }
    
    score <- colSums(gtau_out$G)
    sum_score <- sum(score[1:2])
    
    inv_EI <- ginv(gtau_out$I)[1:2,1:2]  ## J''
    ell <- rep(1, 2)
    variance <- sqrt(c(t(ell) %*% ginv(inv_EI) %*% ell))
    
    result[i] <- sum_score/variance
    if (i %% 100 == 0) print(i)
    
    if(result[i] > 0.6){
      
      ###############
      ## Bootstrap step
      K = nrow(y)
      grid = 50
      score_b <- c()
      
      for(b in 1:B) {
        boot_id = sample(1:param_grid[num, 3], param_grid[num, 3], replace = TRUE)
        yB = y[boot_id, ]
        x1B <- sample_data$x1[boot_id, ]
        x2B <- sample_data$x2[boot_id, ]
        zB <- sample_data$z[boot_id, ]
        
        logLB <- function(paramB) {
          beta1B <- paramB[c(1, 2)]
          beta2B <- paramB[c(3, 4)]
          gammaB <- paramB[c(5, 6)]
          wB <- paramB[7]
          
          d <- 1 - exp(-1)
          
          mu1B <- c(exp(x1B %*% beta1B)); mu2B <- c(exp(x2B %*% beta2B))
          c1B <- exp(mu1B * exp(-1) - mu1B); c2B <- exp(mu2B * exp(-1) - mu2B)
          
          phiB <- c(exp(zB %*% gammaB) / (1 + exp(zB %*% gammaB)))
          indB <- (yB[, 1] == 0 & yB[, 2] == 0)
          
          v1 = exp(zB %*% gammaB) + dpois(x = 0, lambda = mu1B) * dpois(x = 0, lambda = mu2B) * (1 + wB * (1 - c1B) * (1 - c2B))
          v1 <- ifelse(v1 <= 0, 1e-15, v1)
          v1 <- sum(log(v1) * indB)
          
          v21 <- (dpois(x = yB[, 1], lambda = mu1B, log = T) + dpois(x = yB[, 2], lambda = mu2B, log = T)) * !indB
          v22 <- c((1 + wB * (exp(-yB[, 1]) - c1B) * (exp(-yB[, 2]) - c2B)) * !indB)
          v22 <- ifelse(v22 <= 0, log(1e-15), log(v22)) * !indB
          v2 <- (v21 + v22)*!indB
          logbp <- sum(v1) + sum(v2) - sum(log(1+exp(zB %*% gammaB)))
          
          return(-logbp)
        }
        
        paramB <- try(optim(par = rep(0, 7), fn = logLB, method = "L-BFGS-B")$par, silent = TRUE)
        
        if (class(paramB) == "try-error") {
          next
        } else {
          
        }
        
        beta_iB <- paramB[c(1, 3)]
        beta_cB <- paramB[c(2, 4)]
        gamB <- paramB[(c(5, 6))]
        wB <- paramB[7]
        
        
        ## Parameter
        muB <- cbind(exp(x1B %*% c(beta_iB[1], beta_cB[1])), exp(x2B %*% c(beta_iB[2], beta_cB[2])))
        
        ## Gradient
        gtau_outB = tau_gradient(yB, cbind(x1B[, 2], x2B[, 2]), zB[, 2], muB, gamB, wB)
        if (class(ginv(gtau_outB$I))[[1]] == "try-error") {
          next
        } else {
          
        }
        scoreB <- colSums(gtau_outB$G)
        
        # scoreB[1:2] <- ifelse(scoreB[1:2] < 0, 0, scoreB[1:2])
        
        sum_scoreB <- sum(scoreB[1:2])
        
        inv_EIB <- ginv(gtau_outB$I)[1:2,1:2]  ## J''
        ellB <- rep(1, 2)
        varianceB <- sqrt(c(t(ellB) %*% ginv(inv_EIB) %*% ellB))
        
        score_b[b] <- (sum_scoreB-sum_score)/varianceB
        
      }
      
      B_ctic <- quantile(score_b, c(0.99, .95, .9, 0.8), na.rm = TRUE)
      boot_cric[i, ] <- ifelse(result[i] > B_ctic, 1, 0)
    } else{
      boot_cric[i, ] <- c(0, 0)
    }
    
  }
  out = list(result = result, boot_cric = boot_cric)
# }
  
  return(out)
}

result_list3 <- parLapply(cl, X = 1:40, please)

stopCluster(cl)


## 결과  
n_grid <- c(200, 500)
g2_grid <- c(-1.994, -0.373, 0.705, 1.589, 2.4)
w_grid <- c(-1, 0, 1, 2)

param_grid <- expand.grid(w_grid, g2_grid, n_grid)

qval = c(0.05, 0.01)
upp = qchisq(qval, 2, lower.tail = F)/4 + qchisq(qval, 1, lower.tail = F)/2
upp = qchisq(qval, 2, lower.tail = F)
result_list[[1]]$result

out = list()
# 스코어 결과
for(ii in 1:40){
  cat("\n n =", param_grid[ii, 3], "gam2 =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  cat("estimated prob =", round(colMeans(result_list3[[ii]]$boot_cric, na.rm = T), 3), '\n')
  out[[ii]] = c(param_grid[ii, 3], param_grid[ii, 2], param_grid[ii, 1], round(colMeans(result_list3[[ii]]$boot_cric, na.rm = T), 3))
}

# library(xtable)
xtable(do.call("rbind", out))


# 붓스트랩 결과
for(ii in 1:45){
  cat("\n n =", param_grid[ii, 3], "gam2 =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  aa = rowSums(sapply(1:1000, function(j) result_list[[ii]][[2]][j,] < result_list[[ii]][[1]][j]), na.rm = T)/1000
  print(aa)
  cat("\n")
}

save(result_list, file = "~/tau_score.RData")
