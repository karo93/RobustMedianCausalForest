# BLB-style CATE estimation using honest causal trees ----
# Based on Athey & Wager (2016, 2018) and GRF (2019)

blb_cate <- function(
        Y,
        W,
        W.hat,
        X,
        Xnew = NULL,
        split.Rule = split.R,
        cv.option = cv.opt,
        blb_options = .inference$blb,
        ...) {

    ## TRY
    #Y = estim_dataset$Y[idx != i]
    #W = W[idx != i]
    #W.hat = W.hat[idx != i]
    #X = X[idx != i, ]
    #Xnew = X[idx == i, ]
    #split.Rule = split.R
    #cv.option = cv.opt
    # blb_options = .inference$blb

    # -----

    # retrieving option parameters
    B = blb_options$ntrees.bag
    ell = blb_options$nbag
    s_ratio = blb_options$s_ratio
    minsize = blb_options$minsize_blb

    if (is.null(Xnew)) Xnew <- X

    n      <- nrow(X)
    d      <- ncol(X)
    n_new  <- nrow(Xnew)
    n_half <- floor(n / 2)

    # Number of little bags (groups)
    G <- floor(B / ell)
    if (G < 2) stop("Need at least two little bags")

    # Storage:
    # For each group g, we store tree-level predictions at Xnew
    # Dimension: n_new × ell
    group_preds <- vector("list", G)

    # Step 1: grow trees inside half-sample little bags ----
    for (g in seq_len(G)) {
        # --- draw half-sample Hg (theoretical H_g) ---
        Hg <- sample.int(n, n_half, replace = FALSE)

        Yg <- Y[Hg]
        Wg <- W[Hg]
        Xg <- X[Hg, , drop = FALSE]

        # storage for tree-level predictions
        preds_g <- matrix(NA, n_new, ell)

        for (b in seq_len(ell)) {

            # --- draw subsample I_b \in H_g ---
            # smaller subsample from half-sample
            s <- floor(n_half*s_ratio)
            idx <- sample.int(n_half, s, replace = FALSE)

            # treatment balance check (avoid C-level crashes)
            if (length(unique(Wg[idx])) < 2){
                stop("Treatment-control unbalanced.")
            }

            train_df <- as.data.frame(cbind(Y=Yg[idx], Xg[idx, , drop = FALSE]))

            # with try Catch if causalForrest runs into problems
            preds_g[, b] <- tryCatch( {
                # train one tree on subsample s drawn from half-sample H_g
                forest_half <- robustCausalTree::causalForest(
                    as.formula(paste("Y  ~", paste(names(train_df[, -1]), collapse = "+"), sep="")),
                    data = train_df,
                    treatment = Wg[idx], 	# W.res not used here as it just takes W \in {0,1} for now and no continuous W.
                    double.Sample = TRUE, # honest trees
                    split.Rule = split.Rule,
                    cv.option = cv.option,
                    num.trees = 1,
                    nodesize= minsize,# A target for the minimum number of observations in each tree leaf.
                    ncolx = d, # total number of covariates
                    ncov_sample= min(round(d/2), 5) # Number of covariates randomly sampled to build each tree in the forest. At max 5 variables are considered to speed up training.
                    )

                # predict for Xnew
                pr <- predict(forest_half, newdata = Xnew, predict.all = TRUE)

                # save individual estimands
                pr$individual},

                error =  function(e) {
                    # print(glue::glue('g = {g}; and b = {b}: Forest did not run')) # print warning
                    NA                           # return value - must be last!
                }

            ) # And of tryCatch() function

        }

        group_preds[[g]] <- preds_g
    }

    # Step 2: ANOVA decomposition (Section 4.1 in GRF paper) ----

    # Group means:
    group_means <- sapply(
        group_preds,
        function(M) rowMeans(M, na.rm = TRUE)
    )
    # dimension: n_new × G

    # Overall mean:
    overall_mean <- rowMeans(group_means, na.rm = TRUE)

    # ---- Between-group variance (half-sampling variance) ----
    between_var <- rowMeans(
        (group_means - overall_mean)^2,
        na.rm = TRUE
    )

    # ---- Within-group variance (Monte Carlo noise) ----
    within_var <- rowMeans(
        sapply(group_preds, function(M) {
            rowMeans((M - rowMeans(M, na.rm = TRUE))^2, na.rm = TRUE)
        }),
        na.rm = TRUE
    )

    # ---- H-BLB variance estimator (ANOVA formula)  ----
    H_hat_raw <- ((G / (G - 1))  * between_var) - within_var / (ell - 1) # mean should be close to zero, otherwise increase B
    # Bayesian ANOVA shrinkage (Bayesian ANOVA formula)
    # shrinks more when G is small
    # converges to the classical estimator as B approaches infinity
    # H_hat <-  between_var) - within_var / ell)
    H_hat <- pmax(H_hat_raw, 0)   # Truncate at zero (posterior support)

    # Step 3: V-BLB estimator: pseudo-outcome for curvature ----

    # pseudo-outcome for curvature
    curv_y <- (W - W.hat)^2

    # honest regression forest
    rf_V <- grf::regression_forest(X, curv_y, honesty = TRUE)

    V_hat <- predict(rf_V, Xnew)$predictions
    V_inv <- 1/V_hat

    # Step 4: Final variance estimate (see Eq (16) in GRF paper) ----

    Sigma_hat <- V_inv^2 * H_hat
    se_hat    <- sqrt(Sigma_hat)
    uncor_se_hat    <- sqrt(H_hat)


    list(
        tau_se         = se_hat, # corrected SE
        tau_se_uncor   = uncor_se_hat, # uncorrected SE
        H_hat          = H_hat,
        between_var    = between_var, # ist nicht korrigiert
        within_var     = within_var,
        V_hat          = V_hat #, mean_H_hat_raw = mean(H_hat_raw)
    )

}
