rm(list = ls())

library(dplyr)

health <- readxl::read_xlsx("C:/Users/jieun/Dropbox/BZINB/zero-zero/app_data/Aus_health/health.xlsx")

Y <- health[c("DVISITS", "NONDOCCO")] %>% 
  as.matrix()
X <- cbind(1, health[2:12]) %>% 
  as.matrix() %>%
  as.numeric %>%
  matrix(., nrow = 5190, 12)

logH0 <- function(param) {
  beta1 <- param[1:12]
  beta2 <- param[13:24]
  gamma <- param[25:36]
  w <- param[37]
  
  d <- 1 - exp(-1)
  
  ind <- (Y[, 1] == 0 & Y[, 2] == 0)
  
  mu1 <- exp(X %*% beta1)
  mu2 <- exp(X %*% beta2)
  
  phi <- exp(X %*% gamma) / (1 + exp(X %*% gamma))
  
  c1 <- exp(-d * mu1);
  c2 <- exp(-d * mu2)

  V1 <- exp(X %*% gamma) + dpois(0, lambda = mu1) * dpois(0, lambda = mu2) * (1 + w * (1 - c1) * (1 - c2))  
  V1 <- ifelse(V1 <= 0, 1e-15, V1)
  v1 <- log(V1) * ind
  v21 <- (dpois(x = Y[, 1], lambda = mu1, log = T) + dpois(x = Y[, 2], lambda = mu2, log = T)) * !ind
  v22 <- (1 + w * (exp(-Y[, 1]) - c1) * (exp(-Y[, 2]) - c2)) 
  v22 <- ifelse(v22 < 1e-10, 1e-10, v22)
  v22 <- log(v22) * !ind
  v2 <- v21+v22
  
  logbp <- sum(v1) +  sum(v2) - sum(log(1 + exp(X %*% gamma)))
  return(-logbp)
}

logH1 <- function(param) {
  beta1 <- param[1:12]
  beta2 <- param[13:24]
  gamma <- param[25:36]
  w <- param[37]
  tau <- param[38:39]
  
  d <- 1 - exp(-1)
  
  ind <- (Y[, 1] == 0 & Y[, 2] == 0)
  
  mu1 <- exp(X %*% beta1)
  mu2 <- exp(X %*% beta2)
  
  phi <- exp(X %*% gamma) / (1 + exp(X %*% gamma))
  
  c1 <- (1 + d * mu1 * tau[1])^{-1/tau[1]}
  c2 <- (1 + d * mu2 * tau[2])^{-1/tau[2]}
  
  v1 <- log(phi + (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]} * (1 + w * (1 - c1) * (1 - c2)))[ind]
  v2 <- log((1 - phi) * dnbinom(x = Y[, 1], mu = mu1, size = 1/tau[1]) * dnbinom(x = Y[, 2], mu = mu2, size = 1/tau[2]) * (1 + w * (exp(-Y[, 1]) - c1) * (exp(-Y[, 2]) - c2)))[!ind]
  logbp <- sum(v1) +  sum(v2)
  return(-logbp)
}


# initial parameter estimation.
beta1_init <- glm(Y[, 1] ~ X[, 2:12], family = "poisson")$coef
beta2_init <- glm(Y[, 2] ~ X[, 2:12], family = "poisson")$coef
zero_init <- glm(ifelse(rowSums(Y) == 0, 1, 0) ~ X[, 2:12], family = "binomial")$coef

param_init1 <- c(beta1_init, beta2_init, zero_init, .5) %>% 
  as.numeric()
# param_init2 <- c(beta1_init, beta2_init, zero_init, .5, .5, .5) %>% 
#   as.numeric()

## Under H0
param_H0 <- optim(par = param_init1, fn = logH0, method = "L", control = list(maxit = 1e+6), hessian = TRUE)
param_H0$par / round(sqrt(diag(solve(param_H0$hessian))), 3)

## Under H1
# param_H1 <- optim(par = param_init2, fn = logH1, method = "L", control=list(maxit = 1e+6), hessian = TRUE)
# param_H1
## King and Wu test

## Indicator 
ind <- (rowSums(Y) == 0)

d <- 1 - exp(-1)

beta1 <- param_H0$par[1:12]
beta2 <- param_H0$par[13:24]
gamma <- param_H0$par[25:36]
w <- param_H0$par[37]

mu1 <- exp(X %*% beta1)
mu2 <- exp(X %*% beta2)

c1 <- exp(-d * mu1);
c2 <- exp(-d * mu2);

D <- (1 + w * (exp(-Y[, 1]) - c1) * (exp(-Y[, 2]) - c2))

phi <- exp(X %*% gamma) / (1 + exp(X %*% gamma))
psi <- phi / (1 - phi)

grad_mat <- matrix(0, 5190, 39)

## TAU
a1 = sapply(1:nrow(Y), function(i){ 
  sv = 0
  if(Y[i, 1] > 0){
    for(v in 1:Y[i, 1]){
      sv = sv + (Y[i,1] - v)
    }
  }
  return (sv)
})

a2 = sapply(1:nrow(Y), function(i){ 
  sv = 0
  if(Y[i, 2] > 0){
    for(v in 1:Y[i, 2]){
      sv = sv + (Y[i,2] - v)
    }
  }
  return (sv)
})

for (i in 1:12) {
  grad_mat[, i] <- (X[, i] * exp(X %*% gamma) / (exp(X %*% gamma) + exp(-mu1 - mu2) * D)) * ind -
    X[, i] * exp(X %*% gamma) / (1 + exp(X %*% gamma))
  grad_mat[, i + 12] <- (X[, i] * exp(-mu1 - mu2) * mu1 * (w * d * (1 - exp(-d * mu2)) * exp(-d * mu1) - D) / (psi + exp(-mu1 - mu2) * D)) * ind +
    ((Y[, 1] - mu1) * X[, i] + w * d * mu1 * X[, i] * (exp(-Y[, 2]) - exp(-d * mu2)) * exp(-d * mu1) / D) * !ind
  grad_mat[, i + 24] <- (X[, i] * exp(-mu1 - mu2) * mu2 * (w * d * (1 - exp(-d * mu1)) * exp(-d * mu2) - D) / (psi + exp(-mu1 - mu2) * D)) * ind +
    ((Y[, 2] - mu2) * X[, i] + w * d * mu2 * X[, i] * (exp(-Y[, 1]) - exp(-d * mu1)) * exp(-d * mu2) / D) * !ind
}
grad_mat[, 37] <- (exp(-mu1 - mu2) * (1 - exp(-d * mu1)) * (1 - exp(-d * mu2)) / (psi + exp(-mu1 - mu2) * D)) * ind +
  ((exp(-Y[, 1]) - exp(-d * mu1)) * (exp(-Y[, 2]) - exp(-d * mu2)) / D) * !ind
grad_mat[, 38] <- dtau1 <- (exp(-mu1 -mu2) * D * (mu1^2 / 2) / (exp(X %*% gamma) + exp(-mu1 - mu2) * D)) * ind -
  (w * exp(-mu1 - mu2) * exp(-d * mu1) * (1 - exp(-d * mu2)) * (mu1^2 * d / 2)  / (exp(X %*% gamma) + exp(-mu1 - mu2) * D)) * ind +
  (a1 - mu1 * Y[, 1] + mu1^2/ 2 - w * (exp(-Y[, 2]) - c2) * exp(-d * mu1) * (d^2 * mu1^2 / 2) / (1 + D)) * !ind
grad_mat[, 39] <- dtau2 <- (exp(-mu1 -mu2) * D * (mu2^2 / 2) / (exp(X %*% gamma) + exp(-mu1 - mu2) * D)) * ind -
  (w * exp(-mu1 - mu2) * exp(-d * mu2) * (1 - exp(-d * mu1)) * (mu2^2 * d / 2)  / (exp(X %*% gamma) + exp(-mu1 - mu2) * D)) * ind +
  (a2 - mu2 * Y[, 2] + mu2^2/ 2 - w * (exp(-Y[, 1]) - c1) * exp(-d * mu2) * (d^2 * mu2^2 / 2) / (1 + D)) * !ind

King_score <- colSums(grad_mat)

I <- matrix(0, 39, 39)
for (i in 1:39) {
  I[i, ] <- colSums(grad_mat * grad_mat[, i])
}

J <- MASS::ginv(I)[38:39, 38:39]
ell <- rep(1, 2)
JJ <- sqrt(c(t(ell) %*% MASS::ginv(J) %*% ell))
sum_score <- sum(King_score[38:39])
King <- sum_score / JJ

## LR test
# LR_test <- 2 * (param_H0$value - param_H1$value)

## Boot

# i = 1
B = 1000
result <- c()
score_b <- vector(length = B)

library(parallel)

numCores <- detectCores() - 1
cl <- makeCluster(numCores)



boot_ft <- function(num) {
  library(dplyr)
  set.seed(num)
  
  health <- readxl::read_xlsx("C:/Users/jieun/Dropbox/BZINB/zero-zero/app_data/Aus_health/health.xlsx")
  
  Y <- health[c("DVISITS", "NONDOCCO")] %>% 
    as.matrix()
  X <- cbind(1, health[2:12]) %>% 
    as.matrix() %>%
    as.numeric %>%
    matrix(., nrow = 5190, 12)
  
  sum_score <- 1743.261
  
  boot_id <- sample(1:nrow(health), nrow(health), replace = TRUE)
  
  YB = Y[boot_id, ]
  XB <- X[boot_id, ]
  
  logLB <- function(paramB) {
    beta1B <- paramB[1:12]
    beta2B <- paramB[13:24]
    gammaB <- paramB[25:36]
    wB <- paramB[37]
    
    d <- 1 - exp(-1)
    
    indB <- (YB[, 1] == 0 & YB[, 2] == 0)
    
    mu1B <- exp(XB %*% beta1B)
    mu2B <- exp(XB %*% beta2B)
    
    phiB <- exp(XB %*% gammaB) / (1 + exp(XB %*% gammaB))
    
    c1B <- exp(-d * mu1B);
    c2B <- exp(-d * mu2B)
    
    V1 <- exp(XB %*% gammaB) + dpois(0, lambda = mu1B) * dpois(0, lambda = mu2B) * (1 + wB * (1 - c1B) * (1 - c2B))
    V1 <- ifelse(V1 <= 0, 1e-15, V1)
    v1 <- log(V1) * indB
    v21 <- (dpois(x = YB[, 1], lambda = mu1B, log = T) + dpois(x = YB[, 2], lambda = mu2B, log = T)) * !indB
    v22 <- (1 + wB * (exp(-YB[, 1]) - c1B) * (exp(-YB[, 2]) - c2B)) 
    v22 <- ifelse(v22 <= 0, 1e-10, v22)
    v22 <- log(v22) * !indB
    v2 <- v21+v22
    V3 <- 1 + exp(XB %*% gammaB)
    V3 <- ifelse(V3 <= 0, 1e-10, V3)
    v3 <- log(V3)

    logbp <- sum(v1) +  sum(v2) - sum(V3)
    return(-logbp)
  }
  
  beta1_initB <- glm(YB[, 1] ~ XB[, 2:12], family = "poisson")$coef
  beta2_initB <- glm(YB[, 2] ~ XB[, 2:12], family = "poisson")$coef
  zero_initB <- glm(ifelse(rowSums(YB) == 0, 1, 0) ~ XB[, 2:12], family = "binomial")$coef
  
  param_init1B <- c(beta1_initB, beta2_initB, zero_initB, .5) %>% 
    as.numeric()
  
  ## Under H0
  param_H0B <- optim(par = param_init1B, fn = logLB, method = "L", control=list(maxit = 1e+6))
  
  indB <- (rowSums(YB) == 0)
  
  d <- 1 - exp(-1)
  
  beta1B <- param_H0B$par[1:12]
  beta2B <- param_H0B$par[13:24]
  gammaB <- param_H0B$par[25:36]
  wB <- param_H0B$par[37]
  
  mu1B <- exp(XB %*% beta1B)
  mu2B <- exp(XB %*% beta2B)
  
  c1B <- exp(-d * mu1B);
  c2B <- exp(-d * mu2B);
  
  DB <- (1 + wB * (exp(-YB[, 1]) - c1B) * (exp(-YB[, 2]) - c2B))
  
  phiB <- exp(XB %*% gammaB) / (1 + exp(XB %*% gammaB))
  psiB <- phiB / (1 - phiB)
  
  
  grad_matB <- matrix(0, 5190, 39)
  
  ## TAU
  a1B = sapply(1:nrow(YB), function(i){ 
    sv = 0
    if(YB[i, 1] > 0){
      for(v in 1:YB[i, 1]){
        sv = sv + (YB[i,1] - v)
      }
    }
    return (sv)
  })
  
  a2B = sapply(1:nrow(YB), function(i){ 
    sv = 0
    if(YB[i, 2] > 0){
      for(v in 1:YB[i, 2]){
        sv = sv + (YB[i,2] - v)
      }
    }
    return (sv)
  })
  
  for (i in 1:12) {
    grad_matB[, i] <- (XB[, i] * exp(XB %*% gammaB) / (exp(XB %*% gammaB) + exp(-mu1B - mu2B) * DB)) * indB -
      XB[, i] * exp(XB %*% gammaB) / (1 + exp(XB %*% gammaB))
    grad_matB[, i + 12] <- (XB[, i] * exp(-mu1B - mu2B) * mu1B * (wB * d * (1 - exp(-d * mu2B)) * exp(-d * mu1B) - DB) / (psiB + exp(-mu1B - mu2B) * DB)) * indB +
      ((YB[, 1] - mu1B) * XB[, i] + wB * d * mu1B * XB[, i] * (exp(-YB[, 2]) - exp(-d * mu2B)) * exp(-d * mu1B) / DB) * !indB
    grad_matB[, i + 24] <- (XB[, i] * exp(-mu1B - mu2B) * mu2B * (wB * d * (1 - exp(-d * mu1B)) * exp(-d * mu2B) - DB) / (psiB + exp(-mu1B - mu2B) * DB)) * indB +
      ((YB[, 2] - mu2B) * XB[, i] + wB * d * mu2B * XB[, i] * (exp(-YB[, 1]) - exp(-d * mu1B)) * exp(-d * mu2B) / DB) * !indB
  }
  grad_matB[, 37] <- (exp(-mu1B - mu2B) * (1 - exp(-d * mu1B)) * (1 - exp(-d * mu2B)) / (psiB + exp(-mu1B - mu2B) * DB)) * indB +
    ((exp(-YB[, 1]) - exp(-d * mu1B)) * (exp(-YB[, 2]) - exp(-d * mu2B)) / DB) * !indB
  grad_matB[, 38] <- dtau1B <- (exp(-mu1B -mu2B) * DB * (mu1B^2 / 2) / (exp(XB %*% gammaB) + exp(-mu1B - mu2B) * DB)) * indB -
    (wB * exp(-mu1B - mu2B) * exp(-d * mu1B) * (1 - exp(-d * mu2B)) * (mu1B^2 * d / 2)  / (exp(XB %*% gammaB) + exp(-mu1B - mu2B) * DB)) * indB +
    (a1B - mu1B * YB[, 1] + mu1B^2/ 2 - wB * (exp(-YB[, 2]) - c2B) * exp(-d * mu1B) * (d^2 * mu1B^2 / 2) / (1 + DB)) * !indB
  grad_matB[, 39] <- dtau2B <- (exp(-mu1B -mu2B) * DB * (mu2B^2 / 2) / (exp(XB %*% gammaB) + exp(-mu1B - mu2B) * DB)) * indB -
    (wB * exp(-mu1B - mu2B) * exp(-d * mu2B) * (1 - exp(-d * mu1B)) * (mu2B^2 * d / 2)  / (exp(XB %*% gammaB) + exp(-mu1B - mu2B) * DB)) * indB +
    (a2B - mu2B * YB[, 2] + mu2B^2/ 2 - wB * (exp(-YB[, 1]) - c1B) * exp(-d * mu2B) * (d^2 * mu2B^2 / 2) / (1 + DB)) * !indB
  
  
  IB <- matrix(0, 39, 39)
  for (i in 1:39) {
    IB[i, ] <- colSums(grad_matB * grad_matB[, i])
  }
  
  JB <- MASS::ginv(IB)[38:39, 38:39]
  ellB <- rep(1, 2)
  JJB <- sqrt(c(t(ellB) %*% MASS::ginv(JB) %*% ellB))
  King_scoreB <- colSums(grad_matB)
  sum_scoreB <- sum(King_scoreB[38:39])
  KingB <- (sum_scoreB - sum_score) / JJB
  return(KingB)
}

xx = parLapply(cl = cl, X = 1:1000, fun = boot_ft)
save2 = do.call("c", xx2)

save = c(save, KingB)
hist(save, breaks=100)
