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
#' @param minsize_estim minimum size of leaf during forest training
#' @param split.R splitting rule (one of "CT", "MSD", "MAD", "LMS")
#' @param cv.opt cross-valdiation option (one of "CT", "MSD", "MAD", "LMS")
#' @param ATE_estimation # should ATE be computed and given out via crossfitted AIPW. Default TRUE.
#' @param .inference The variance specific options to use for variance estimation. This
#'   must be the result from a call to [inference_options()].
#' @returns A list with two elements. The first one is \code{pred_cate} which holds
#' the predicted CATE and the corresponding standard error for each unit.
#' The second element is \code{pred_ate} which holds the ATE and its corresponding
#' standard error, calculated based on debiased CATEs through an AIPW estimator.
#'
#' @rdname robustCausalForest
#' @export
#' @aliases robustCausalForest
#'
robustCausalForest_BLB_IJ_sim <- function(Y, # outcome variable
                               W, # treatment variable
                               X, # covariate matrix
                               Y.hat = NULL,
                               W.hat = NULL,
                               Y.res = NULL,
                               W.res = NULL,
                               k = 4, # number of folds for cross-fitting
                               idx, # vector of indices for cross-fitting
                               ntrees =  4e3, # number of trees in a forest
                               minsize_estim = 5, # minimum number of obs in every leaf for training/estimation
                               split.R = "CT",
                               cv.opt = "CT",
                               ATE_estimation = TRUE,
                               .inference = inference_options(), # Options for inference (variance) estimation
                       ...){

  on.exit(gc(), add = TRUE)  # clean on exit

  # tries
  # Y = df$sim_dataset[[1]]$Y
  # W = df$sim_dataset[[1]]$W
  #  X = df$sim_dataset[[1]] %>% dplyr::select(-c(Y, W, tau))
  # Y.hat = df$Y.hat[[1]]
  # W.hat = df$W.hat[[1]]
  # Y.res = df$Y.res[[1]]
  # W.res = df$W.res[[1]]
  # k = 4
  # idx = df$idx[[1]]
  # ntrees = 4e3
  # split.R = "CT"
  # cv.opt = "CT"
  # ATE_estimation = TRUE
  # .inference = df$inference[[1]]
  # minsize_estim = 5

  # start main function ----

  # Sample ids for cross-fitting
  D = ncol(X) # number of covariats in X
  N = nrow(X) # number of obs. in dataset

  # Initialize storing vectors and matrices
  tau.hat <- matrix(rep(NaN, N*ntrees), ncol=ntrees, nrow=N)
  inbag.inf <- vector("list", k)
  # A list to contain IJ and BLB var estimation

  .inference$varianceMethod <- c(.inference$varianceMethod, paste0(.inference$varianceMethod, "_uncor"))

  if (length(.inference$varianceMethod) > 0) {
    tau.se.hat <- setNames(vector("list", length(.inference$varianceMethod)),
                          .inference$varianceMethod)
    for (method in .inference$varianceMethod) {
    tau.se.hat[[method]] <- rep(NaN, N)
    }
  }

  # Prepare dataset for causalForest training ----

  # residualized or unresidualized outcome
  if(!is.null(Y.res)){
    # use Y-residuals to center out nuisance
    estim_dataset <- as.data.frame(cbind(Y=Y.res, X))
  } else if(is.null(Y.res)) {
    estim_dataset <- as.data.frame(cbind(Y=Y, X))
    warning("Y.res. not specified. Estimation of CATEs proceeds uncentered using Y.")
  } else (
    warning("Specify Y.res. Default is Y.res=NULL and estimation proceeds uncentered.")
  )

  # For CATE; center TRUE/FALSE -> center res
  # BLB oder ATE -> center hats

  # I have this in the sim outside the robustCausalForest function oder beides? !!!
  # estimate W.hat, Y.hat, if needed (for BLB and ATE via AIPW)
  if("BLB" %in% .inference$varianceMethod || ATE_estimation == TRUE){
    if(is.null(W.hat) || is.null(Y.hat) ){

      Y.hat <- rep(NaN, length(Y))
      W.hat <- rep(NaN, length(W))


      center_data <- data.frame(W=W, X)

      for (i in 1:k){
        # zu GRF umschreiben
        W.model <- SimDesign::quiet(randomForest::randomForest(W~., center_data[idx != i, ], ntree=2000))
        W.hat[idx == i] <- SimDesign::quiet(predict(W.model, newdata = X[idx == i, ]))

        #! what is with res?
      }

      center_data <- cbind(Y=Y, X)

      for (i in 1:k){
        Y.model <- randomForest::randomForest(Y~., center_data[idx != i, ], ntree=2000)
        Y.hat[idx == i] <- predict(Y.model, newdata = X[idx == i, ])
      }

      center_data <- data.frame(Y=Y, W=W, X)

      # new with regression_forest
      for (i in 1:k){
        Y.model <- grf::regression_forest(X = X[idx != i, ], Y = center_data$Y[idx != i], num.trees = ntrees)
        Y.hat[idx == i] <- predict(Y.model, newdata = X[idx == i, ])$predictions
        # Y.res[idx == i] <- center_dataset_Y$Y[idx == i] - Y.hat[idx == i]

        W.model <- grf::regression_forest(X = X[idx != i, ], Y = center_data$W[idx != i], num.trees = ntrees)
        W.hat[idx == i] <- predict(W.model, newdata = X[idx == i, ])$predictions
        # W.res[idx == i] <- center_dataset_W$W[idx == i] - W.hat[idx == i]

      }

    }

  }

  # OOB CATE estimation ----

  for(i in 1:k){
    # Train honest causal forests
    forest_honest <- robustCausalTree::causalForest(
      as.formula(paste("Y ~", paste(names(estim_dataset[, -1]), collapse = "+"), sep="")),
      data = estim_dataset[idx != i, ],
      treatment = W[idx != i], 	# W.res not used here as it just takes W \in {0,1} for now and no continuous W.
      double.Sample = TRUE, # honest trees
      split.Rule = split.R,
      cv.option = cv.opt,
      num.trees = ntrees,
      nodesize= minsize_estim,# A target for the minimum number of observations in each tree leaf.
      ncolx = D, # total number of covariates
      ncov_sample= min(round(D/2), 5) # Number of covariates randomly sampled to build each tree in the forest. At max 5 variables are considered.
    )

    # ── Robinson correction: estimate V̂(x) on training fold ──────────────────
    curv_y <- (W[idx != i] - W.hat[idx != i])^2 # squared residuals on propensity score

    rf_V  <- grf::regression_forest(
      X[idx != i, ],
      curv_y,
      honesty = TRUE
    )

    V_hat <- predict(rf_V, X[idx == i, ])$predictions
    V_inv <- 1 / V_hat

    # V_inv <- 1 / pmax(V_hat, 1e-6)  # maybe implement
    # ─────────────────────────────────────────────────────────────────────────


    # OOB predictions
    preds <- predict(forest_honest, newdata =  X[idx == i, ], predict.all = TRUE)
    inbag.inf[[i]] <- forest_honest$inbag
    tau.hat[idx == i, ] <- preds$individual

    # Bootstrap inference function ----
    if ("IJ" %in% .inference$varianceMethod) {
      # Estimate variance via IJ
      # https://rdrr.io/cran/ranger/src/R/infinitesimalJackknife.R
      H_hat_IJ <- ranger:::rInfJack(
        pred = preds$individual,
        inbag = forest_honest$inbag,
        calibrate = TRUE)$var.hat

      # Step 2: Robinson correction (Athey, Tibshirani & Wager 2019, Eq. 16)
      # Var[τ̂(x)] = V̂(x)^{-2} * Ĥ(x)

      tau.se.hat$IJ_uncor[idx == i] <- sqrt(H_hat_IJ)
      tau.se.hat$IJ[idx == i] <- sqrt(V_inv^2 * H_hat_IJ)
    } else{tau.se.hat$IJ[idx == i] <- NULL}

    if ("BLB" %in% .inference$varianceMethod) {
      blb_se <- blb_cate(Y = estim_dataset$Y[idx != i], # Why Y from here and X not?
                         W = W[idx != i],
                         W.hat = W.hat[idx != i],
                         X = X[idx != i, ],
                         Xnew = X[idx == i, ],
                         split.Rule = split.R,
                         cv.option = cv.opt,
                         blb_options = .inference$blb)

      tau.se.hat$BLB[idx == i] <- blb_se$tau_se
      tau.se.hat$BLB_uncor[idx == i] <- blb_se$tau_se_uncor
    } else{tau.se.hat$BLB[idx == i] <- NULL}

    # print(glue::glue("Finished run {i} of {k} OOB folds (k), at: {format(Sys.time(), '%d.%m.%Y %H:%M:%S')}."))

    rm(forest_honest, preds); gc()
  }

  # prepare savings ----

  # old
  #cate_est_old <- data.frame(estimate=rowMeans(tau.hat),
  #                       # New for different methods
  #                       std.err= sqrt(tau.var.hat$IJ)) %>%
  #  mutate(id = row_number())

  # new with tibble Maybe change to more vars and not nested list
  cate_est <- tibble::as_tibble(tau.se.hat) %>%
    dplyr::rename_with(~glue::glue("std.err_{.}")) %>%
    dplyr::mutate(estimate = rowMeans(tau.hat)) %>%
    mutate(id = row_number())

  ### ATE estimation ----

  if(ATE_estimation){

    # list for storing
    ate_est <- list()
    # Preparation for grf:::average_treatment_effect()
    forest_obj <- list()
    forest_obj$Y.orig <- as.vector(Y)
    forest_obj$W.orig <- as.vector(W)
    forest_obj$W.hat <- W.hat
    forest_obj$Y.hat <- Y.hat
    # it is not used in avg_treat_eff, is it?
    forest_obj$tau.hat.pointwise <- cate_est$estimate

    ate_est <- avg_treat_eff(forest=forest_obj,
                             target.sample = "all",
                             method = "AIPW")


    } else if(ATE_estimation==FALSE){

    ate_est <- NULL

  } else {
    # warning to the top?
    warning("ATE_estimation needs to be TRUE/FALSE.")
  }

  rm(tau.hat)

  return(list(pred_cate = cate_est, pred_ate=ate_est))
}


# Average Treatment Effect function
# Simplified version of grf::average_treatment_effect


avg_treat_eff <- function(forest,
                          target.sample = "all",
                          method = "AIPW",
                          # not sure here
                          subset = NULL,
                          debiasing.weights = NULL,
                          compliance.score = NULL,
                          num.trees.for.weights = 500) {

  # try
  # forest <- forest_obj
  # print('Inside avg_treat_eff function')

  clusters <- 1:NROW(forest$Y.orig)
  observation.weight <- obs_weights(forest)
# I do
  subset <- valid_subset(forest, subset)
  subset.clusters <- clusters[subset]
  subset.weights <- observation.weight[subset]

  if (length(unique(subset.clusters)) <= 1) {
    stop("The specified subset must contain units from more than one cluster.")
  }

  # What happens with NULL in debiasing.weights?
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


    DR.scores <- compute_DR_scores(forest,
                               subset = subset,
                               # where does the debiasing.weigths coming from? should be NULL right?
                               debiasing.weights = NULL, # debiasing.weights,
                               compliance.score = NULL, # compliance.score,
                               num.trees.for.weights = num.trees.for.weights)


    tau.hat <- weighted.mean(DR.scores, subset.weights)
    # where is tau.hat coming from? here only the forest Object is present
    sigma2.hat <- .sigma2.hat(DR.scores, tau.hat, subset.clusters, subset.weights)
    return(c(estimate = tau.hat, std.err = sqrt(sigma2.hat)))
  }
}

# doesnt the subset.cluster and weigths are missing?
.sigma2.hat <- function(DR.scores, tau.hat, subset.clusters, subset.weights, ..) {
  correction.clust <-
    Matrix::sparse.model.matrix(~ factor(subset.clusters) + 0,
                                transpose = TRUE) %*%
    (sweep(as.matrix(DR.scores), 2, tau.hat, "-") * subset.weights)

  # effective number of samples(clusters) is all with > 0 weight.
  n.adj <- sum(rowsum(subset.weights, subset.clusters) > 0)

  Matrix::colSums(correction.clust^2) / sum(subset.weights)^2 * n.adj / (n.adj - 1)
}

obs_weights <- function(forest) {
  # Case: No given sample.weights

  raw.weights <- rep(1, NROW(forest$Y.orig))

  return (raw.weights / sum(raw.weights))
}

valid_subset <- function(forest, subset) {
  if (is.null(subset)) {
    # if subset is null creat a vector
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
  return(subset)
}


compute_DR_scores <- function(forest,
                          subset = NULL,
                          debiasing.weights = NULL,
                          num.trees.for.weights = 500,
                          ...) {

  subset <- valid_subset(forest, subset)
  W.orig <- forest$W.orig[subset]
  W.hat <- forest$W.hat[subset]
  Y.orig <- forest$Y.orig[subset]
  Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- forest$tau.hat.pointwise

  binary.W <- all(forest$W.orig %in% c(0, 1))

  if (is.null(debiasing.weights)) {
    if (binary.W) {
      debiasing.weights <- (W.orig - W.hat) / (W.hat * (1 - W.hat))
    }# else {
      # Start by learning debiasing weights if needed.
      # The goal is to estimate the variance of W given X. For binary treatments,
      # we get a good implicit estimator V.hat = e.hat (1 - e.hat), and
      # so this step is not needed. Note that if we use the present CAPE estimator
      # with a binary treatment and set V.hat = e.hat (1 - e.hat), then we recover
      # exactly the AIPW estimator of the CATE.
    #   clusters <- if (length(forest$clusters) > 0) {
    #     forest$clusters
    #   } else {
    #     1:length(forest$Y.orig)
    #   }
    #   variance_forest <- grf::regression_forest(forest$X.orig,
    #                                             (forest$W.orig - forest$W.hat)^2,
    #                                             clusters = clusters,
    #                                             sample.weights = forest$sample.weights,
    #                                             num.trees = num.trees.for.weights,
    #                                             ci.group.size = 1,
    #                                             seed = forest$seed,
    #                                             num.threads = forest$num.threads)
    #   V.hat <- predict(variance_forest)$predictions
    #   debiasing.weights.all <- (forest$W.orig - forest$W.hat) / V.hat
    #   debiasing.weights <- debiasing.weights.all[subset]
    # }
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


