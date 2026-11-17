gen_data <- function(sim = 1, n = 1000, ntest = 1000, d = 10,
                     indep_covariate, homosc_err, rho) {
  
  sim_id <- rep(1:sim, each = n)
  errdist <- rnorm
  
  X <- Xfun(n * sim, d, indep_covariate, rho)
  colnames(X) <- paste("X", seq(d), sep = "")
  
  Y0 <- rep(0, n * sim)
  tau <- taufun(X)
  std <- sdfun(X, homosc_err)
  Y1 <- tau + std * errdist(n * sim)
  
  ps <- psfun(X)
  T_ind <- as.numeric(runif(n * sim) < ps)
  
  Yobs <- Y0
  Yobs[T_ind == 1] <- Y1[T_ind == 1]
  
  Xtest <- Xfun(ntest * sim, d, indep_covariate, rho)
  colnames(Xtest) <- paste("Xtest", seq(d), sep = "")
  
  Y0test <- rep(0, ntest * sim)
  tautest <- taufun(Xtest)
  stdtest <- sdfun(Xtest, homosc_err)
  Y1test <- tautest + stdtest * errdist(ntest * sim)
  pstest <- psfun(Xtest)
  
  dat <- tibble(
    sim_id,
    Y1,
    Y0,
    tau,
    ps,
    T_ind,
    Yobs,
    as.tibble(X),
    Y1test,
    Y0test,
    tautest,
    pstest,
    as.tibble(Xtest)
  )
  
  arrange(dat, sim_id)
  
  return(dat)
}


# Wager and Athey (2018), equation (27)

Xfun <- function(n, d, indep_covariate, rho) {
  
  if (indep_covariate == TRUE) {
    Xfunmat <- matrix(runif(n * d), nrow = n, ncol = d)
    return(Xfunmat)
    
  } else if (indep_covariate == FALSE) {
    X <- matrix(rnorm(n * d), nrow = n, ncol = d)
    fac <- rnorm(n)
    X <- X * sqrt(1 - rho) + fac * sqrt(rho)
    Xfunmat <- pnorm(X)
    return(Xfunmat)
    
  } else {
    stop("Neither indep covariates nor corr covariates chosen.")
  }
}

taufun <- function(X) {
  2 / (1 + exp(-12 * (X[, 1] - 0.5))) *
    2 / (1 + exp(-12 * (X[, 2] - 0.5)))
}

sdfun <- function(X, homosc_err) {
  
  if (homosc_err == TRUE) {
    return(rep(1, nrow(X)))
    
  } else if (homosc_err == FALSE) {
    return(-log(X[, 1] + 1e-9))
    
  } else {
    stop("Neither homosc. nor heterosc. errors chosen.")
  }
}

psfun <- function(X) {
  (1 + pbeta(X[, 1], 2, 4)) / 4
}


user_Dev <- function(Y, X, Xtest, T_ind) {
  
  X <- as.data.frame(X)
  p <- ncol(X)
  Xtest <- as.data.frame(Xtest)
  
  data_YXT <- data.frame(Y = Y, X, T_ind)
  names(Xtest) <- names(data_YXT)[c(-1, -ncol(data_YXT))]
  
  f <- paste(
    names(data_YXT)[1], "~",
    paste(names(Xtest), collapse = " + ")
  )
  
  fit <- causalTree::causalForest(
    as.formula(f),
    data = data_YXT,
    treatment = data_YXT$T_ind,
    split.Rule = "user",
    cv.option = "user",
    split.Honest = TRUE,
    cv.Honest = TRUE,
    split.Bucket = FALSE,
    bucketNum = 5,
    split.alpha = 0.5,
    nodesize = 5,
    ncov_sample = floor(p / 3),
    ncolx = p,
    minsize = 2,
    propensity = 0.5,
    num.trees = 1000
  )
  
  as.numeric(predict(fit, Xtest))
}


CT_Dev <- function(Y, X, Xtest, T_ind) {
  
  X <- as.data.frame(X)
  p <- ncol(X)
  Xtest <- as.data.frame(Xtest)
  
  data_YXT <- data.frame(Y = Y, X, T_ind)
  names(Xtest) <- names(data_YXT)[c(-1, -ncol(data_YXT))]
  
  f <- paste(
    names(data_YXT)[1], "~",
    paste(names(Xtest), collapse = " + ")
  )
  
  fit <- causalTree::causalForest(
    as.formula(f),
    data = data_YXT,
    treatment = data_YXT$T_ind,
    split.Rule = "CT",
    cv.option = "CT",
    split.Honest = TRUE,
    cv.Honest = TRUE,
    split.Bucket = FALSE,
    bucketNum = 5,
    split.alpha = 0.5,
    nodesize = 5,
    ncov_sample = floor(p / 3),
    ncolx = p,
    minsize = 2,
    propensity = 0.5,
    num.trees = 1000
  )
  
  as.numeric(predict(fit, Xtest))
}

