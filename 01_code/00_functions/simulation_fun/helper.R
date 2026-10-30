## DGP of Athey and Wager (2018) for their Procedure 1 ----

# -------------------------------
# Data-generating process
# -------------------------------

sigma <- function(x, sig_variant = 2, ...) {
    if(sig_variant == 1) return(1 + 1 / (1 + exp(-20 * (x - 1/3))))
    if(sig_variant == 2) return(2 / (1 + exp(-12 * (x - 1/2))))
    if(!sig_variant %in% 1:2) stop('variant not defined')
}

tau_function <- function(X, tau_variant, sig_variant = 2, ...) {
  if(tau_variant == 'complex_linear'){
    n_col <- NCOL(X)
    beta_0 <- runif(n_col, 1, 30)
    beta_1 <- runif(n_col, 1, 30)
    tau <- X %*% (beta_1 - beta_0) # transposed?
    return(as.numeric(tau))
  }
  if(tau_variant == 'complex_non_linear'){
    return(sigma(X[, 1], sig_variant = sig_variant) *
           sigma(X[, 2], sig_variant = sig_variant))
  }
  if(tau_variant == 'cnl_interaction'){
    tau <-  6 * sin(2 * X[, 1]) +
      3 * (X[, 2] + 3) * X[, 3] +
      9 * tanh(0.5 * X[, 4]) +
      3 * X[, 5] * (2 * as.numeric(X[, 4] > 0) - 1) +
      3 * X[, 6] + 2 * X[, 7] + X[, 8] - 2 * X[, 9] - 4 * X[, 10]
    return(as.numeric(tau))
  }
  if(tau_variant == 'no_effect'){
    return(rep(0, NROW(X)))
  }
  if(!tau_variant %in% c('complex_linear', 'complex_non_linear', 'no_effect', 'cnl_interaction')) {
    stop('Tau variant not definied! must be one of: complex_linear, complex_non_linear, or no_effect')
  }
}

heavy_noise <- function(n, df = 3, ...) {
  rt(n, df = df) / sqrt(df / (df - 2))
}

laplace_noise <- function(n, noise_parameters = NULL, ...) {
  # default standard laplace
  if(is.null(noise_parameters)){
    loc <- 0
    scale <- 1
  }
  if(sum(names(noise_parameters) %in% c('loc', 'scale')) == 2) {
    loc <- noise_parameters['loc'][[1]]
    scale <- noise_parameters['scale'][[1]]
  }
  if(sum(names(noise_parameters) %in% c('loc', 'scale')) != 2 && !is.null(noise_parameters)){
    stop('Noise distirbution parameters must be provided in named list!')
  }

  return(rLaplace(n, mu = loc, b = scale))
}


cauchy_noise <- function(n, noise_parameters = NULL, ...){
  # default standard cauchy
  if(is.null(noise_parameters)){
    loc <- 0
    scale <- 1
  }
  if(sum(names(noise_parameters) %in% c('loc', 'scale')) == 2) {
    loc <- noise_parameters['loc'][[1]]
    scale <- noise_parameters['scale'][[1]]
  }
  if(sum(names(noise_parameters) %in% c('loc', 'scale')) != 2 && !is.null(noise_parameters)){
    stop('Noise distirbution parameters must be provided in named list!')
  }

    return(rcauchy(n, location = loc, scale = scale))
}


small_region_smooth <- function(X, r = 1.2) {
  as.integer(X[,1]^2 + X[,2]^2 > r^2)
}

simulate_data <- function(n, d,
                          noise, # one of c("normal", "heavy", "cauchy", "laplace")
                          noise_parameters = NULL, # list of parameters
                          tau_effect, # one of c("het", "het.sparse")
                          tau_variant, # one of c('complex_linear', 'complex_non_linear', 'no_effect')
                          sigma_variant = 1, # 1 = 1 + 1 / (1 + exp(-20 * (x - 1/3))) # 2 =
                          x_dist = 'uniform', # either uniform or normal
                          ...) {
  if(x_dist == 'uniform') X <- matrix(runif(n * d), n, d)   # all runif
  if(x_dist == 'normal') X <- matrix(rnorm(n * d), n, d)   # all normal
  if(!x_dist %in% c('normal', 'uniform')) stop('Distribution of `X` must be either normal or uniform. Default: uniform.')
  W <- rbinom(n, 1, 0.5)            # randomized treatment

  # tau
  if (tau_effect=="het"){
    tau <- tau_function(X, tau_variant, sigma_variant)
  } else if(tau_effect=="het.sparse"){
    # Sparse extreme responders
    magnitude <- 10
    region <- small_region_smooth(X)
    tau <- tau_function(X, tau_variant, sigma_variant) + magnitude * region
  } else{
    stop("Choose right 'noise' option.")
  }

  # obs.
  if (noise=="normal"){
    Y0 <- rnorm(n) # m(x) = 0, just noise independent of x
    Y1 <- tau + Y0
    Y <- W * Y1 + (1-W)*Y0
  } else if(noise=="heavy"){
    Y0 <- heavy_noise(n)
    Y1 <- tau + Y0 # m(x) = 0, just noise independent of x
    Y <- W * Y1 + (1-W)*Y0
  }else if(noise=="cauchy"){
      Y0 <- cauchy_noise(n, noise_parameters = noise_parameters)
      Y1 <- tau + Y0 # m(x) = 0, just noise independent of x
      Y <- W * Y1 + (1 - W)*Y0
  }else if(noise=="laplace"){
    Y0 <- laplace_noise(n, noise_parameters = noise_parameters)
    Y1 <- tau + Y0 # m(x) = 0, just noise independent of x
    Y <- W * Y1 + (1 - W)*Y0
  } else{
    stop("Choose right 'noise' option.")
  }

  return(data.frame(X = X, W = W, Y = Y, tau = tau))
}

gen_dataset <- function(num_dataset, N, D, tau_effect, tau_variant, sigma_variant = NULL, noise, noise_parameters, x_dist = 'uniform', ...){

  # set seeds
  seeds <- sample.int(1e6, num_dataset)

  # parse noise parameter for naming
  if(is.null(noise_parameters)){
    noise_parameter_p <- 'sd.'}
  else {
    noise_parameter_p <- paste0(stringr::str_sub(names(noise_parameters), 1, 1), noise_parameters, collapse = "")
  }

  tau_variant_p <- ifelse(tau_variant == "complex_linear", 'cl',
                          ifelse(tau_variant == "complex_non_linear", 'cnl', 'ne')
                          ) #c('complex_linear', 'complex_non_linear', 'no_effect')

  sigma_variant_p <- ifelse(is.null(sigma_variant), "", paste0('_', sigma_variant))


  # generate folder to saves files
  path <- here::here('00_sim_data', tau_variant, glue::glue('{tau_effect}{sigma_variant_p}'), glue('{noise}_{noise_parameter_p}'), glue::glue('N.{N}'), x_dist)
  dir.create(file.path(path), recursive = TRUE, showWarnings = FALSE)

  for(i in seq_len(num_dataset)){

    # get seed
    seed <- seeds[i]

    # generate data
    sim_dataset <- simulate_data(n=N, d=D, noise=noise, tau_effect=tau_effect,
                                 tau_variant = tau_variant,
                                 sigma_variant = sigma_variant,
                                 seed = seed, x_dist = x_dist)

    saveRDS(sim_dataset, file = glue::glue('{path}/ds{formatC(i, width = 3, flag = "0")}.N{N/1000}.tau{tau_variant_p}_{tau_effect}{sigma_variant_p}.X{stringr::str_sub(x_dist, 1, 1)}.noise{noise}_{noise_parameter_p}.rds'))

    if(i %% 10 == 0) print(paste("Datset ", i, " of ", num_dataset, " finished."))
  }
}

####################################
####--- evaluation functions ---####
####################################

rmse_eval <- function(pred, actual){

  sqrt(mean((pred-actual)^2))
}

absbias_eval <- function(pred, actual){

  abs(pred-actual)

}

coverage_width_eval <- function(estim, std.err, true_tau){
  # 95% interval
  low <- estim - std.err*1.96
  up <- estim + std.err*1.96
  coverage <- ifelse(up >= true_tau & low <= true_tau, 1, 0)
  width <- up - low

  return(
    list(coverage = coverage,
         width = width)
    )
}


