

simulate_correlation_matrix <- function(dim, alpha) {
  betaparam <- 1 / alpha
  P <- matrix(nrow = dim, ncol = dim)
  S <- diag(dim)
  
  for (k in 1:(dim - 1)) {
    for (i in (k + 1):dim) {
      P[k, i] <- (rbeta(1, betaparam, betaparam) - 0.5) * 2
      p <- P[k, i]
      
      if (k > 1) {
        for (l in (k - 1):1) {
          p <- p * sqrt((1 - P[l, i]^2) * (1 - P[l, k]^2)) +
            P[l, i] * P[l, k]
        }
      }
      
      S[k, i] <- p
      S[i, k] <- p
    }
  }
  
  permutation <- sample(seq_len(dim))
  S[permutation, permutation]
}

Xfunc_2 <- function(d, alph, n) {
  if (alph == 0) {
    X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  } else {
    corr.mat <- simulate_correlation_matrix(dim = d, alpha = alph)
    X <- mvtnorm::rmvnorm(n, sigma = corr.mat)
  }
  
  colnames(X) <- paste0("X", seq_len(d))
  X
}

sigma_fun <- function(x) {
  2 / (1 + exp(-12 * (x - 0.5)))
}

gen_data_1 <- function(n = 1000, d = 10, d_mu = 10,
                       prop.score = 0.5, alph = 0) {
  X <- Xfunc_2(d = d, alph = alph, n = n)
  
  beta_0 <- runif(d, min = 1, max = 10)
  beta_1 <- runif(d, min = 1, max = 10)
  
  mu_0 <- as.vector(X[, 1:d_mu, drop = FALSE] %*% beta_0[1:d_mu])
  mu_1 <- as.vector(X[, 1:d_mu, drop = FALSE] %*% beta_1[1:d_mu])
  CATE_true <- mu_1 - mu_0
  
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size = 1, prob = prop.score)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  tibble(
    Y1 = Y1,
    Y0 = Y0,
    Yobs = Yobs,
    T_ind = T_ind,
    CATE_true = CATE_true,
    ps = prop.score
  ) %>%
    bind_cols(as_tibble(X))
}

gen_data_2 <- function(n = 1000, d = 10,
                       prop.score = 0.5, alph = 0) {
  X <- Xfunc_2(d = d, alph = alph, n = n)
  
  mu_0 <- -0.5 * sigma_fun(X[, 1]) * sigma_fun(X[, 2])
  mu_1 <- 0.5 * sigma_fun(X[, 1]) * sigma_fun(X[, 2])
  CATE_true <- mu_1 - mu_0
  
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size = 1, prob = prop.score)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  tibble(
    Y1 = Y1,
    Y0 = Y0,
    Yobs = Yobs,
    T_ind = T_ind,
    CATE_true = CATE_true,
    ps = prop.score
  ) %>%
    bind_cols(as_tibble(X))
}

gen_data_3 <- function(n = 1000, d = 10, d_mu = 10,
                       prop.score = 0.5, alph = 0) {
  X <- Xfunc_2(d = d, alph = alph, n = n)
  
  beta <- runif(d, min = 1, max = 10)
  
  mu_0 <- as.vector(X[, 1:d_mu, drop = FALSE] %*% beta[1:d_mu])
  mu_1 <- mu_0
  CATE_true <- mu_1 - mu_0
  
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size = 1, prob = prop.score)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  tibble(
    Y1 = Y1,
    Y0 = Y0,
    Yobs = Yobs,
    T_ind = T_ind,
    CATE_true = CATE_true,
    ps = prop.score
  ) %>%
    bind_cols(as_tibble(X))
}

gen_data_4 <- function(n = 1000, d = 10, alph = 0) {
  X <- Xfunc_2(d = d, alph = alph, n = n)
  
  beta_0 <- function(X) {
    100 + 4 * X[, 1] + X[, 2] - 3 * X[, 3]
  }
  
  CATE_fun <- function(X) {
    6 * sin(2 * X[, 1]) +
      3 * (X[, 2] + 3) * X[, 3] +
      9 * tanh(0.5 * X[, 4]) +
      3 * X[, 5] * (2 * as.numeric(X[, 4] > 0) - 1) +
      3 * X[, 6] + 2 * X[, 7] + X[, 8] -
      2 * X[, 9] - 4 * X[, 10]
  }
  
  CATE_true <- CATE_fun(X)
  mu_0 <- beta_0(X) - 0.5 * CATE_true
  mu_1 <- beta_0(X) + 0.5 * CATE_true
  
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size = 1, prob = 0.5)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  tibble(
    Y1 = Y1,
    Y0 = Y0,
    Yobs = Yobs,
    T_ind = T_ind,
    CATE_true = CATE_true,
    ps = 0.5
  ) %>%
    bind_cols(as_tibble(X))
}

sim <- 1000
n <- 1000
d <- 10

scenario <- "S4"

generators <- list(
  S1 = gen_data_1,
  S2 = gen_data_2,
  S3 = gen_data_3,
  S4 = gen_data_4
)

generator <- generators[[scenario]]

output_dir <- here::here("00_data", scenario, "00_Clean")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1234)

train <- replicate(
  sim,
  generator(n = n, d = d, alph = 0),
  simplify = FALSE
)

train <- bind_rows(train, .id = "sim_id")
train$sim_id <- as.integer(train$sim_id)

saveRDS(
  train,
  file.path(output_dir, paste0(scenario, "_train.rds"))
)

test <- replicate(
  sim,
  generator(n = n, d = d, alph = 0),
  simplify = FALSE
)

test <- bind_rows(test, .id = "sim_id")
test$sim_id <- as.integer(test$sim_id)

saveRDS(
  test,
  file.path(output_dir, paste0(scenario, "_test.rds"))
)

