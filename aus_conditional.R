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
  tau2 <- param[38]

  d <- 1 - exp(-1)

  mu1 <- c(exp(X %*% beta1)); mu2 <- c(exp(X %*% beta2))
  c1 <- exp(-d * mu1);
  c2 <- (1 + d * mu2 * tau2)^{-tau2^{-1}}

  ind <- (Y[, 1] == 0 & Y[, 2] == 0)

  v1 <- log(exp(X %*% gamma) + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))[ind]
  v2 <- (dpois(x = Y[, 1], lambda = mu1, log = T) + dnbinom(x = Y[, 2], mu = mu2, size = 1/tau2, log = T))[!ind]
  v3 <- (1 + w * (exp(-Y[, 1]) - c1) * (exp(-Y[, 2]) - c2))

  v3 <- ifelse(v3 <= 0, log(1e-10), log(v3))[!ind]

  logbp <- sum(v1) + sum(v2) + sum(v3) - sum(log(1+exp(X %*% gamma)))

  return(-logbp)
}

# initial parameter estimation.
beta1_init <- glm(Y[, 1] ~ X[, 2:12], family = "poisson")$coef
beta2_init <- glm(Y[, 2] ~ X[, 2:12], family = "poisson")$coef
zero_init <- glm(ifelse(rowSums(Y) == 0, 1, 0) ~ X[, 2:12], family = "binomial")$coef

param_init1 <- c(beta1_init, beta2_init, zero_init, .5, .5) %>%
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
tau2 <- param_H0$par[38]

mu1 <- exp(X %*% beta1)
mu2 <- exp(X %*% beta2)

c1 <- exp(-d * mu1);
c2 <- (1 + d * mu2 * tau2)^{-tau2^{-1}}

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
  grad_mat[, i] <- (((exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * mu1 * X[, i] * (w * d * (1 - c2) * c1 - (1 + w * (1 - c1) * (1 - c2)))) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind) +
    (((Y[, 1] - mu1) * X[, i] + w * d * mu1 * X[, i] * (exp(-Y[, 2]) - c2) * c1 / (1 + w * (exp(-Y[, 1]) - c1) * (exp(-Y[, 2]) - c2))) * !ind)
  grad_mat[, i + 12] <- ((exp(-mu1) * mu2 * X[, i] * (1 + tau2 * mu2)^{-1/tau2} * (w * d * (1 - c1) * (1 + d * tau2 * mu2)^{-1/tau2 - 1} - (1 + tau2 * mu2)^{-1} * (1 + w * (1 - c1) * (1 - c2))) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 -c2)))) * ind) +
    (((Y[, 2] - mu2) * X[, i] / (1 + tau2 * mu2) + w * d * mu2 * X[, i] * (exp(-Y[, 1]) - c1) * (1 + d * tau2 * mu2)^{-1/tau2 - 1} / (1 + w * (exp(-Y[, 1])- c1) * (exp(-Y[, 2]) - c2))) * !ind)
  grad_mat[, i + 24] <- (psi *  X[, i] / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind - psi * X[, i] /(1 + psi)
}
grad_mat[, 37] <- ((exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)) * (tau2^{-2} * log(1 + tau2 * mu2) - mu2 * tau2^{-1} / (1 + tau2 * mu2)) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind) -
  ((w * (1 - c1) * exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * c2 * (tau2^{-2} * log(1 + d * tau2 * mu2) - d * mu2 * tau2^{-1} / (1 + d * tau2 * mu2)) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind) +
  ((a2 - Y[, 2] * mu2 / (1 + tau2 * mu2) + tau2^{-2} * log(1 + tau2 * mu2) - tau2^{-1} * mu2 / (1 + tau2 * mu2)) * !ind) -
  ((w * (exp(-Y[, 1]) - c1) * c2 * (tau2^{-2} * log(1 + d * tau2 * mu2) - d * mu2 * tau2^{-1} / (1 + d * tau2 * mu2)) / (1 + w * (exp(-Y[, 1]) - c1) * (exp(-Y[, 2]) - c2))) * !ind)
grad_mat[, 38] <- (exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 - c1) * (1 - c2) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind +
  ((exp(-Y[, 1]) - c1) * (exp(-Y[, 2]) - c2) / (1 + w * (exp(-Y[, 1]) - c1) * (exp(-Y[, 2]) - c2))) * !ind
grad_mat[, 39] <- (((exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)) * mu1^2/2) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind)  -
  (((w * c1 * (1 - c2) * exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * d^2 * mu1^2 / 2) / (psi + exp(-mu1) * (1 + tau2 * mu2)^{-1/tau2} * (1 + w * (1 - c1) * (1 - c2)))) * ind) +
  ((a1 - Y[, 1] * mu1 + mu1^2 / 2 - w * c1 * (exp(-Y[, 2]) - c2) * (d * mu1)^2 / 2 / (1 + w * (exp(-Y[, 1]) - c1) * (exp(-Y[, 2]) - c2))) * !ind)

score <- colSums(grad_mat)

I <- matrix(0, 39, 39)
for (i in 1:39) {
  I[i, ] <- colSums(grad_mat * grad_mat[, i])
}

I11 <- I[1:38, 1:38]
I22 <- diag(I)[39]
I12 <- I[1:38, 39]

result <- score[39] / sqrt(I22 - t(I12) %*% MASS::ginv(I11) %*% (I12))
## LR test
# LR_test <- 2 * (param_H0$value - param_H1$value)

## Boot

# i = 1
B = 1000
result <- c()
score_b <- vector(length = B)

library(parallel)

cl <- makeCluster(60)

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

  logLB <- function(param) {
    beta1B <- param[1:12]
    beta2B <- param[13:24]
    gammaB <- param[25:36]
    wB <- param[37]
    tau2B <- param[38]

    d <- 1 - exp(-1)

    mu1B <- exp(XB %*% beta1B)
    mu2B <- exp(XB %*% beta2B)
    phiB <- exp(XB %*% gammaB) / (1 + exp(XB %*% gammaB))

    c1B <- exp(-d * mu1B);
    c2B <- (1 + d * mu2B * tau2B)^{-tau2B^{-1}}

    indB <- (YB[, 1] == 0 & YB[, 2] == 0)

    v1 <- exp(XB %*% gammaB) + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B))
    v1 <- log(v1)[indB]
    v2 <- (dpois(x = YB[, 1], lambda = mu1B, log = T) + dnbinom(x = YB[, 2], mu = mu2B, size = 1/tau2B, log = T))[!indB]
    v3 <- (1 + wB * (exp(-YB[, 1]) - c1B) * (exp(-YB[, 2]) - c2B))
    v3 <- ifelse(v3 <= 0, log(1e-10), log(v3))[!indB]
    v4 <- 1 + exp(XB %*% gammaB)
    v4 <- ifelse(v4 <= 0, log(1e-10), log(v4))

    logbp <- sum(v1) + sum(v2) + sum(v3) - sum(v4)
    return(-logbp)
  }

  beta1_initB <- glm(YB[, 1] ~ XB[, 2:12], family = "poisson")$coef
  beta2_initB <- glm(YB[, 2] ~ XB[, 2:12], family = "poisson")$coef
  zero_initB <- glm(ifelse(rowSums(YB) == 0, 1, 0) ~ XB[, 2:12], family = "binomial")$coef

  param_init1B <- c(beta1_initB, beta2_initB, zero_initB, .5, .5) %>%
    as.numeric()

  ## Under H0
  param_H0B <- try({optim(par = param_init1B, fn = logLB, method = "L",
                          lower = c(rep(-10, 37), 1e-10), upper = c(rep(10, 37), 10),
                          control=list(maxit = 200))})

  if(class(param_H0B) == "try-error"){
    KingB <- NA
  } else{
  indB <- (rowSums(YB) == 0)

  d <- 1 - exp(-1)

  beta1B <- param_H0B$par[1:12]
  beta2B <- param_H0B$par[13:24]
  gammaB <- param_H0B$par[25:36]
  wB <- param_H0B$par[37]
  tau2B <- param_H0B$par[38]

  mu1B <- exp(XB %*% beta1B)
  mu2B <- exp(XB %*% beta2B)

  c1B <- exp(-d * mu1B);
  c2B <- (1 + d * mu2B * tau2B)^{-tau2B^{-1}}

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
    grad_matB[, i] <- (((exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * mu1B * XB[,i] * (wB * d * (1 - c2B) * c1B - (1 + wB * (1 - c1B) * (1 - c2B)))) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB) +
      (((YB[, 1] - mu1B) * XB[,i] + wB * d * mu1B * XB[,i] * (exp(-YB[, 2]) - c2B) * c1B / (1 + wB * (exp(-YB[, 1]) - c1B) * (exp(-YB[, 2]) - c2B))) * !indB)
    grad_matB[, i + 12] <- ((exp(-mu1B) * mu2B * XB[,i] * (1 + tau2B * mu2B)^{-1/tau2B} * (wB * d * (1 - c1B) * (1 + d * tau2B * mu2B)^{-1/tau2B - 1} - (1 + tau2B * mu2B)^{-1} * (1 + wB * (1 - c1B) * (1 - c2B))) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 -c2B)))) * indB) +
      (((YB[, 2] - mu2B) * XB[,i] / (1 + tau2B * mu2B) + wB * d * mu2B * XB[,i] * (exp(-YB[, 1]) - c1B) * (1 + d * tau2B * mu2B)^{-1/tau2B - 1} / (1 + wB * (exp(-YB[, 1])- c1B) * (exp(-YB[, 2]) - c2B))) * !indB)
    grad_matB[, i + 24] <- (psiB *  XB[,i] / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB - psiB * XB[,i] /(1 + psiB)
  }
  grad_matB[, 37] <- ((exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)) * (tau2B^{-2} * log(1 + tau2B * mu2B) - mu2B * tau2B^{-1} / (1 + tau2B * mu2B)) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB) -
    ((wB * (1 - c1B) * exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * c2B * (tau2B^{-2} * log(1 + d * tau2B * mu2B) - d * mu2B * tau2B^{-1} / (1 + d * tau2B * mu2B)) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB) +
    ((a2B - YB[, 2] * mu2B / (1 + tau2B * mu2B) + tau2B^{-2} * log(1 + tau2B * mu2B) - tau2B^{-1} * mu2B / (1 + tau2B * mu2B)) * !indB) -
    ((wB * (exp(-YB[, 1]) - c1B) * c2B * (tau2B^{-2} * log(1 + d * tau2B * mu2B) - d * mu2B * tau2B^{-1} / (1 + d * tau2B * mu2B)) / (1 + wB * (exp(-YB[, 1]) - c1B) * (exp(-YB[, 2]) - c2B))) * !indB)
  grad_matB[, 38] <- (exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 - c1B) * (1 - c2B) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB +
    ((exp(-YB[, 1]) - c1B) * (exp(-YB[, 2]) - c2B) / (1 + wB * (exp(-YB[, 1]) - c1B) * (exp(-YB[, 2]) - c2B))) * !indB
  grad_matB[, 39] <- (((exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)) * mu1B^2/2) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB)  -
    (((wB * c1B * (1 - c2B) * exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * d^2 * mu1B^2 / 2) / (psiB + exp(-mu1B) * (1 + tau2B * mu2B)^{-1/tau2B} * (1 + wB * (1 - c1B) * (1 - c2B)))) * indB) +
    ((a1B - YB[, 1] * mu1B + mu1B^2 / 2 - wB * c1B * (exp(-YB[, 2]) - c2B) * (d * mu1B)^2 / 2 / (1 + wB * (exp(-YB[, 1]) - c1B) * (exp(-YB[, 2]) - c2B))) * !indB)

  scoreB <- colSums(grad_matB, na.rm = T)

  IB <- matrix(0, 39, 39)
  for (i in 1:39) {
    IB[i, ] <- colSums(grad_matB * grad_matB[, i])
  }

  I11B <- IB[1:38, 1:38]
  I22B <- diag(IB)[39]
  I12B <- IB[1:38, 39]

  boot <- (scoreB[39] - 450.7114) / sqrt(I22B - t(I12B) %*% MASS::ginv(I11B) %*% (I12B))

  }

  return(boot)
}

xx = parLapply(cl = cl, X = 1:1000, fun = boot_ft)
save = do.call("c", xx)
hist(save, breaks = 100)
abline(v = King, col = "blue")
quantile(save, 0.95, na.rm = T)

save = c(save, KingB)
hist(save, breaks=100)
