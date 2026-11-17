# data generation künzel et al. 
gen_data_1 = function(n=1000, # number of observations
                      d=20, # number of independent variables
                      d_mu=10, # number of independent variables, which effect mu, default=8
                      prop.score=0.5, 
                      alph = 0# if 0, X's are independent, if alph>0, X's are correlated  
){
  
  
  X <- Xfunc_2(d=d, alph=alph, n=n)
  beta <- runif(d, min=-5, max=5)
  mu_0 <- X[,1:d_mu] %*% beta[1:d_mu] + as.vector(ifelse(X[, 1] > 0.5, 5, 0))
  CATE_true <- as.vector(ifelse(X[, 2] > 0.1, 8, 0))
  mu_1 <- mu_0 + CATE_true
  
  
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size=1, prob=prop.score)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  dat = tibble(
    Y1=Y1[,1],
    Y0=Y0[,1],
    Yobs,
    T_ind, 
    CATE_true,
    ps=prop.score, 
    as_tibble(X)
  )
  
  return(dat)
  
}


# uses the C-vine method for
#'   simulating correlation matrices.
#'   Lewandowski D, Kurowicka D, Joe H (2009) Generating random correlation matrices based on vines and extended onion method.
#'   https://github.com/soerenkuenzel/causalToolbox/blob/master/R/ExmpleSetups_causaleffects.R
simulate_correlation_matrix <- function(dim, alpha) {
  
  betaparam <- 1 / alpha
  
  P <- matrix(nrow = dim, ncol = dim)           # storing partial correlations
  S <- diag(dim)
  
  for (k in 1:(dim - 1)) {
    for (i in (k + 1):dim) {
      P[k, i] <- rbeta(1, betaparam, betaparam) # sampling from beta
      P[k, i] <- (P[k, i] - 0.5) * 2     # linearly shifting to [-1, 1]
      p <- P[k, i]
      if (k > 1) {
        for (l in (k - 1):1) {
          # converting partial correlation to raw correlation
          p <- p * sqrt((1 - P[l, i] ^ 2) * (1 - P[l, k] ^ 2)) + P[l, i] *
            P[l, k]
          # p  = p * sqrt((1-P(l,i)^2)*(1-P(l,k)^2)) + P(l,i)*P(l,k);
        }
      }
      S[k, i] <- p
      S[i, k] <- p
    }
  }
  
  # permuting the variables to make the distribution permutation-invariant
  permutation <- sample(1:dim)
  S <- S[permutation, permutation]
  return(S)
}

Xfunc_2 <- function(d, alpha, n){
  corr.mat <- simulate_correlation_matrix(dim=d, alpha=alph)
  x <- mvtnorm::rmvnorm(n, sigma=corr.mat)
  colnames(x) <- paste("X", seq(d), sep="")
  return(x)
}

#######################
gen_data_2 = function(n=1000, # number of observations
                      d=20, # number of independent variables
                      d_mu=10, # number of independent variables, which effect mu
                      prop.score=0.5, 
                      alph = 0# if 0, X's are independent, if alph>0, X's are correlated  
){
  
  X <- Xfunc_2(d=d, alph=alph, n=n)
  beta_1 <- runif(d, min=1, max=10)
  beta_0 <- runif(d, min=1, max=10)
  
  mu_0 <- X[,1:d_mu] %*% beta_0[1:d_mu] 
  mu_1 <- X[,1:d_mu] %*% beta_1[1:d_mu] 
  CATE_true <- mu_1 - mu_0 #should be right?
  
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size=1, prob=prop.score)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  
  dat = tibble(
    Y1=Y1[,1],
    Y0=Y0[,1],
    Yobs,
    T_ind, 
    CATE_true=CATE_true[,1],
    ps=prop.score, 
    as_tibble(X)
  )
  
  return(dat)
  
}

gen_data_3 = function(n=1000, # number of observations
                      d=20, # number of independent variables
                      prop.score=0.5, 
                      alph = 0# if 0, X's are independent, if alph>0, X's are correlated  
){
  
  
  X <- Xfunc_2(d=d, alph=alph, n=n)
  
  mu_0 <- -1/2 * sigma_fun(X[, 1]) * sigma_fun(X[, 2])
  mu_1 <- 1/2 * sigma_fun(X[, 1]) * sigma_fun(X[, 2])
  CATE_true <- mu_1 - mu_0 #should be right?
  
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size=1, prob=prop.score)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  
  dat = tibble(
    Y1=Y1,
    Y0=Y0,
    Yobs,
    T_ind, 
    CATE_true,
    ps=prop.score, 
    as_tibble(X)
  )
  
  return(dat)
  
}
sigma_fun <- function(x) {
  y<- 2/((1 + exp(-12*(x - 1/2))))
  return(y)
}

gen_data_4 = function(n=1000, # number of observations
                      d=20, # number of independent variables
                      d_mu=10, # number of independent variables, which effect mu
                      prop.score=0.5, 
                      alph = 0# if 0, X's are independent, if alph>0, X's are correlated  
){
  
  X <- Xfunc_2(d=d, alph=alph, n=n)
  beta <- runif(d, min=1, max=10)
  
  mu_0 <- X[,1:d_mu] %*% beta[1:d_mu] 
  mu_1 <- mu_0
  CATE_true <- mu_1 - mu_0 #should be right?
  
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size=1, prob=prop.score)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  
  dat = tibble(
    Y1=Y1[,1],
    Y0=Y0[,1],
    Yobs,
    T_ind, 
    CATE_true=CATE_true[,1],
    ps=prop.score, 
    as_tibble(X)
  )
  
  return(dat)
  
}
#######################
gen_data_5 = function(n=1000, # number of observations
                      d=100, # number of independent variables
                      alph = 0# if 0, X's are independent, if alph>0, X's are correlated  
){
  
  X <- Xfun(sim, d, indep_covariate, rho)
  colnames(X) <- paste0("X", 1:ncol(X))
  
  #logit_p <- function(x) exp(x) / (1 + exp(x))
  #prop.score <- logit_p(X[, 1] - X[, 2])
  prop.score=0.5
  beta_0 <- function(X) {
    100 + 4 * X[, 1] + X[, 2] - 3 * X[, 3]
  }
  CATE_true <- function(X) {
    6 * sin(2 * X[, 1]) +
      3 * (X[, 2] + 3) * X[, 3] +
      9 * tanh(0.5 * X[, 4]) +
      3 * X[, 5] * (2 * as.numeric(X[, 4] > 0) - 1)+
      3 * X[, 6] + 2 * X[, 7] + X[, 8] - 2 * X[, 9] - 4 * X[, 10]
  }
  mu_0 <- beta_0(X) - 1/2 * CATE_true(X)
  mu_1 <- beta_0(X) + 1/2 * CATE_true(X)
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size=1, prob=prop.score)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  dat = tibble(
    Y1=Y1,
    Y0=Y0,
    Yobs,
    T_ind, 
    CATE_true=CATE_true(X),
    ps=prop.score, 
    as_tibble(X)
  )
  
  return(dat)
  
}

#######################
sim=1000
n=1000

set.seed(1234)
train <- replicate(sim, gen_data_5(n=n, d=10, alph = 0), simplify=FALSE)
sim_id <- tibble(sim_id=rep(1:sim, each=n))
train <- bind_rows(train)
train <- bind_cols(sim_id, train)
saveRDS(object=train, file=here::here("00_data/03_simulation_li_et_al/sim_low/00_Clean/sim_s3_train.rds"))

test <- replicate(sim, gen_data_5(n=n, d=10, alph = 0), simplify=FALSE)
sim_id <- tibble(sim_id=rep(1:sim, each=n))
test <- bind_rows(test)
test <- bind_cols(sim_id, test)
saveRDS(object=test, file=here::here("00_data/03_simulation_li_et_al/sim_low/00_Clean/sim_s3_test.rds"))

######################


