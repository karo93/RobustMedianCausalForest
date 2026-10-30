#' Estimate Heterogeneous Treatment Effect via Causal Forest with Median Splitting Rules
#'
#' Estimate heterogeneous treatment effect based on random forests.
#'
#' @param Y outcome variable
#' @param W treatment variable
#' @param X covariate matrix
#' @param Y.hat cross-fitted estimates of E[Y|X]
#' @param W.hat cross-fitted estimates of E[W|X]
#' @param Y.res outcome residuals
#' @param W.res treatment residuals
#' @param k number of folds for cross-fitting
#' @param idx vector of indices for cross-fitting
#' @param ntrees number of trees in forest
#' @param split.R splitting rule (one of "CT", "MSD", "MAD", "LMS")
#' @param cv.opt cross-valdiation option (one of "CT", "MSD", "MAD", "LMS")
#' @returns A list with two elements. The first one is \code{pred_cate} which holds
#' the predicted CATE and the corresponding standard error for each unit.
#' The second element is \code{pred_ate} which holds the ATE and its corresponding
#' standard error, calculated based on debiased CATEs through an AIPW estimator.
#'
#' @rdname robustCausalForest
#' @export
#' @aliases robustCausalForest
#'
robustCausalForest <- function(Y, # outcome variable
                               W, # treatment variable
                               X, # covariate matrix
                               Y.hat,
                               W.hat,
                               Y.res,
                               W.res,
                               k, # number of folds for cross-fitting
                               idx,
                               ntrees, # number of trees in a forest
                               split.R = "CT",
                               cv.opt = "CT",
                       ...){

  # Sample ids for cross-fitting
  D = ncol(X) # number of covariats in X
  N = nrow(X) # number of obs. in dataset
  #idx <- sample(c(1:k), N, replace = T) # don't sample indices again to do proper cross-fitting

  # Initialize storing vectors and matrices
  tau.hat <- matrix(rep(NaN, N*ntrees), ncol=ntrees, nrow=N)
  inbag.inf <- vector("list", k)
  tau.var.hat <- rep(NaN, N)


  # Prepare dataset for causalForest training ----
  # use Y-residuals to center out nuisance
  estim_dataset <- as.data.frame(cbind(Y.res=Y.res, X))

  # OOB CATE estimation ----

  for(i in 1:k){


    # Train honest causal forests
    forest_honest <- robustCausalTree::causalForest(
      as.formula(paste("Y.res ~", paste(names(estim_dataset[, -1]), collapse = "+"), sep="")),
      data = estim_dataset[idx != i, ],
      treatment = W[idx != i], 	# W.res not used here as it just takes W \in {0,1} for now and no continuous W.
      double.Sample = TRUE, # honest trees
      split.Rule = split.R,
      cv.option = cv.opt,
      num.trees = ntrees,
      nodesize= 5,# A target for the minimum number of observations in each tree leaf.
      ncolx = D, # total number of covariates
      ncov_sample= min(round(D/2), 5) # Number of covariates randomly sampled to build each tree in the forest. At max 5 variables are considered.
    )


    # OOB predictions
    preds <- predict(forest_honest, newdata =  X[idx == i, ], predict.all = TRUE)
    inbag.inf[[i]] <- forest_honest$inbag
    tau.hat[idx == i, ] <- preds$individual
    tau.var.hat[idx == i] <- rInfJack(pred=preds$individual,
                                      forest_honest$inbag,
                                      calibrate = TRUE)$var.hat

    print(paste("finished run ", i, " of k =", k, "OOB folds."))
    print(Sys.time())

  }

  # prepare savings ----
  cate_est <- data.frame(estimate=rowMeans(tau.hat),
                         std.err=sqrt(tau.var.hat)) %>%
    mutate(id = row_number())

  ### ATE estimation ----


  # Create forest-object to estimate doubly-robust ATE from CATEs
  forest_obj <- list()
  forest_obj$Y.orig <- as.vector(Y)
  forest_obj$W.orig <- as.vector(W)
  forest_obj$W.hat <- W.hat
  forest_obj$Y.hat <- Y.hat
  forest_obj$tau.hat.pointwise <- cate_est$estimate
  forest_obj$tau.var.hat.pointwise <- cate_est$std.err

  ate_est <- avg_treat_eff(forest=forest_obj,
                           target.sample = "all",
                           method = "AIPW")



  return(list(pred_cate = cate_est, pred_ate=ate_est))
}


# Average Treatment Effect function

avg_treat_eff <- function(forest,
                          target.sample = "all",
                          method = "AIPW",
                          subset = NULL,
                          debiasing.weights = NULL,
                          compliance.score = NULL,
                          num.trees.for.weights = 500) {


  clusters <- 1:NROW(forest$Y.orig)

  observation.weight <- observation_weights(forest)

  subset <- validate_subset(forest, subset)
  subset.clusters <- clusters[subset]
  subset.weights <- observation.weight[subset]

  if (length(unique(subset.clusters)) <= 1) {
    stop("The specified subset must contain units from more than one cluster.")
  }

  if (!is.null(debiasing.weights)) {
    if (length(debiasing.weights) == NROW(forest$Y.orig)) {
      debiasing.weights <- debiasing.weights[subset]
    } else if (length(debiasing.weights) != length(subset)) {
      stop("If specified, debiasing.weights must be a vector of length n or the subset length.")
    }
  }

  if (method == "AIPW" && target.sample == "all") {
    # This is the most general workflow, that shares codepaths with best linear projection
    # and other average effect estimators.

    .sigma2.hat <- function(DR.scores, tau.hat) {
      correction.clust <- Matrix::sparse.model.matrix(~ factor(subset.clusters) + 0, transpose = TRUE) %*%
        (sweep(as.matrix(DR.scores), 2, tau.hat, "-") * subset.weights)
      n.adj <- sum(rowsum(subset.weights, subset.clusters) > 0) # effective number of samples(clusters) is all with > 0 weight.

      Matrix::colSums(correction.clust^2) / sum(subset.weights)^2 *
        n.adj / (n.adj - 1)
    }


    DR.scores <- get_DR_scores(forest,
                               subset = subset,
                               debiasing.weights = debiasing.weights,
                               compliance.score = compliance.score,
                               num.trees.for.weights = num.trees.for.weights)


    tau.hat <- weighted.mean(DR.scores, subset.weights)
    sigma2.hat <- .sigma2.hat(DR.scores, tau.hat)
    return(c(estimate = tau.hat, std.err = sqrt(sigma2.hat)))
  }
}




observation_weights <- function(forest) {
  # Case: No given sample.weights

  raw.weights <- rep(1, NROW(forest$Y.orig))

  return (raw.weights / sum(raw.weights))
}

validate_subset <- function(forest, subset) {
  if (is.null(subset)) {
    subset <- 1:NROW(forest$Y.orig)
  }
  if (is.logical(subset) && length(subset) == NROW(forest$Y.orig)) {
    subset <- which(subset)
  }
  if (!all(subset %in% 1:NROW(forest$Y.orig))) {
    stop(paste(
      "If specified, subset must be a vector contained in 1:n,",
      "or a boolean vector of length n."
    ))
  }
  subset
}


get_DR_scores <- function(forest,
                          subset = NULL,
                          debiasing.weights = NULL,
                          num.trees.for.weights = 500,
                          ...) {
  subset <- validate_subset(forest, subset)
  W.orig <- forest$W.orig[subset]
  W.hat <- forest$W.hat[subset]
  Y.orig <- forest$Y.orig[subset]
  Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- forest$tau.hat.pointwise

  binary.W <- all(forest$W.orig %in% c(0, 1))

  if (is.null(debiasing.weights)) {
    if (binary.W) {
      debiasing.weights <- (W.orig - W.hat) / (W.hat * (1 - W.hat))
    } else {
      # Start by learning debiasing weights if needed.
      # The goal is to estimate the variance of W given X. For binary treatments,
      # we get a good implicit estimator V.hat = e.hat (1 - e.hat), and
      # so this step is not needed. Note that if we use the present CAPE estimator
      # with a binary treatment and set V.hat = e.hat (1 - e.hat), then we recover
      # exactly the AIPW estimator of the CATE.
      clusters <- if (length(forest$clusters) > 0) {
        forest$clusters
      } else {
        1:length(forest$Y.orig)
      }
      variance_forest <- grf::regression_forest(forest$X.orig,
                                                (forest$W.orig - forest$W.hat)^2,
                                                clusters = clusters,
                                                sample.weights = forest$sample.weights,
                                                num.trees = num.trees.for.weights,
                                                ci.group.size = 1,
                                                seed = forest$seed,
                                                num.threads = forest$num.threads)
      V.hat <- predict(variance_forest)$predictions
      debiasing.weights.all <- (forest$W.orig - forest$W.hat) / V.hat
      debiasing.weights <- debiasing.weights.all[subset]
    }
  } else if (length(debiasing.weights) == length(forest$Y.orig)) {
    debiasing.weights <- debiasing.weights[subset]
  } else if (length(debiasing.weights) != length(subset))  {
    stop("If specified, debiasing.weights must have length n or |subset|.")
  }

  # Form AIPW scores. Note: We are implicitly using the following
  # estimates for the regression surfaces E[Y|X, W=0/1]:
  # Y.hat.0 <- Y.hat - W.hat * tau.hat.pointwise
  # Y.hat.1 <- Y.hat + (1 - W.hat) * tau.hat.pointwise
  Y.residual <- Y.orig - (Y.hat + tau.hat.pointwise * (W.orig - W.hat))

  tau.hat.pointwise + debiasing.weights * Y.residual
}

###############################################################################
# @Descrip - compute CIs for causalForest based on inf. jackknife
# (see randomForestCI and ranger docs for original rInfJack)
###############################################################################


rInfJack = function(pred, inbag, calibrate = TRUE, used.trees = NULL) {
  # compute confidence intervals based on infinitesimalJackknife approach
  # original: https://github.com/swager/randomForestCI/blob/master/R/infinitesimalJackknife.R

  if (is.null(used.trees)) {
    used.trees = 1:ncol(inbag)
  }
  pred = pred[, used.trees, drop=FALSE]

  # Check if sampling without replacement
  no.replacement = (max(inbag) == 1)

  # Extract tree-wise predictions and variable counts from random forest
  B = length(used.trees)
  n = nrow(inbag)
  s = sum(inbag) / ncol(inbag)

  tau.hat = rowMeans(pred)
  pred.centered = pred - rowMeans(pred)

  N = Matrix::Matrix(inbag[, used.trees], sparse = TRUE)
  N.avg = Matrix::rowMeans(N)

  # Compute raw infinitesimal jackknife
  if (as.numeric(B)^2 > as.numeric(n) * as.numeric(nrow(pred))) {

    C = Matrix::tcrossprod(N, pred.centered) -
      Matrix::Matrix(N.avg, nrow(N), 1) %*%
      Matrix::Matrix(rowSums(pred.centered), 1, nrow(pred.centered))
    raw.IJ = Matrix::colSums(C^2) / B^2

  } else {
    # Faster implementation when n is large. Uses the fact that
    # colSums((A - B)^2) = T1 - 2 * T2 + T3,
    # where T1 = diag(A'A), T2 = diag(B'A), and T3 = diag(B'B)

    NTN = Matrix::crossprod(N, N)
    NTNPT_T = Matrix::tcrossprod(pred.centered, NTN)
    T1 = Matrix::rowSums(pred.centered * NTNPT_T)

    RS = rowSums(pred.centered)
    NbarTN = Matrix::crossprod(N.avg, N)
    T2 = RS * Matrix::tcrossprod(NbarTN, pred.centered)

    T3 = sum(N.avg^2) * RS^2
    raw.IJ = as.numeric(T1 - 2 * T2 + T3) / B^2

  }

  # Apply Monte Carlo bias correction
  N.var = mean(Matrix::rowMeans(N^2) - Matrix::rowMeans(N)^2)
  boot.var = rowSums(pred.centered^2) / B
  bias.correction = n * N.var * boot.var / B
  vars = raw.IJ - bias.correction

  # Finite sample correction
  if (no.replacement) {
    variance.inflation = 1 / (1 - mean(inbag))^2
    vars = variance.inflation * vars
  }

  results = data.frame(tau.hat=tau.hat, var.hat=vars)

  if (isTRUE(calibrate) && nrow(results) <= 20) {
    calibrate = FALSE
    warning("Sample size <=20, no calibration performed.")
  }

  # If appropriate, calibrate variance estimates; this step in particular
  # ensures that all variance estimates wil be positive.

  if (calibrate) {
    # Compute variance estimates using half the trees
    calibration.ratio = 2
    n.sample = ceiling(B / calibration.ratio)
    results.ss = rInfJack(pred, inbag, calibrate = FALSE, used.trees = sample(used.trees, n.sample))

    # Use this second set of variance estimates to estimate scale of Monte Carlo noise
    sigma2.ss = mean((results.ss$var.hat - results$var.hat)^2)
    delta = n.sample / B
    sigma2 = (delta^2 + (1 - delta)^2) / (2 * (1 - delta)^2) * sigma2.ss

    # Use Monte Carlo noise scale estimate for empirical Bayes calibration
    results = tryCatch(
      expr = {
        vars.calibrated = calibrateEB(vars, sigma2)
        results$var.hat = vars.calibrated
        results
      },
      error = function(e) {
        warning(sprintf("Calibration failed with error:\n%sFalling back to non-calibrated variance estimates.", e))
        results = rInfJack(pred, inbag, calibrate = FALSE, used.trees = used.trees)
        return(results)
      }
    )
  }

  return(results)
}


# Fit an empirical Bayes prior in the hierarchical model
#     mu ~ G, X ~ N(mu, sigma^2)
#
gfit = function(X, sigma, p = 2, nbin = 1000, unif.fraction = 0.1) {

  xvals = seq(min(min(X) - 2 * sd(X), 0), max(max(X) + 2 * sd(X), sd(X)), length.out = nbin)
  binw = xvals[2] - xvals[1]

  zero.idx = max(which(xvals <= 0))
  noise.kernel = dnorm(xvals / sigma) * binw / sigma

  if (zero.idx > 1) {
    noise.rotate = noise.kernel[c(zero.idx:length(xvals), 1:(zero.idx - 1))]
  } else {
    noise.rotate = noise.kernel
  }

  XX = sapply(1:p, function(j) xvals^j * as.numeric(xvals >= 0))
  neg.loglik = function(eta) {
    g.eta.raw = exp(XX %*% eta) * as.numeric(xvals >= 0)
    if ((sum(g.eta.raw) == Inf) | (sum(g.eta.raw) <= 100 * .Machine$double.eps)) {
      return (1000 * (length(X) + sum(eta^2)))
    }
    g.eta.main = g.eta.raw / sum(g.eta.raw)
    g.eta = (1 - unif.fraction) * g.eta.main +
      unif.fraction * as.numeric(xvals >= 0) / sum(xvals >= 0)
    f.eta = convolve(g.eta, noise.rotate)
    sum(approx(xvals, -log(pmax(f.eta, 0.0000001)), X)$y)
  }

  eta.hat = nlm(neg.loglik, rep(-1, p))$estimate
  g.eta.raw = exp(XX %*% eta.hat) * as.numeric(xvals >= 0)
  g.eta.main = g.eta.raw / sum(g.eta.raw)
  g.eta = (1 - unif.fraction) * g.eta.main +
    unif.fraction * as.numeric(xvals >= 0) / sum(xvals >= 0)

  return(data.frame(x=xvals, g=g.eta))
}

# Bayes posterior estimation with Gaussian noise
#
gbayes = function(x0, g.est, sigma) {
  Kx = dnorm((g.est$x - x0) / sigma)
  post = Kx * g.est$g
  post = post / sum(post)
  sum(post * g.est$x)
}

# Empirical Bayes calibration of noisy variance estimates
#
calibrateEB = function(vars, sigma2) {

  if(sigma2 <= 0 | min(vars) == max(vars)) {
    return(pmax(vars, 0))
  }

  sigma = sqrt(sigma2)
  eb.prior = gfit(vars, sigma)

  if (length(vars >= 200)) {
    # If there are many test points, use interpolation to speed up computations
    calib.x = unique(quantile(vars, q = seq(0, 1, by = 0.02)))
    calib.y = sapply(calib.x, function(xx) gbayes(xx, eb.prior, sigma))
    calib.all = approx(x=calib.x, y=calib.y, xout=vars)$y
  } else {
    calib.all = sapply(vars, function(xx) gbayes(xx, eb.prior, sigma))
  }
  return(calib.all)
}


