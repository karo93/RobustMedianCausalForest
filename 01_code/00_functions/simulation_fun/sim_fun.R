# function
run_sim <- function(dgp_setting, methods = c("CT", "LMS", "MAD", "MSD"),
                    k = 4, ntrees = 4e3, verbose = TRUE, ...){

    on.exit(gc(), add = TRUE)  # clean on exit

    # checking if sim setting already exist ----
    # get noise_parameters
    noise_parameters <- dgp_setting$noise_parameters[[1]]


    # prepare noise parameter for naming
    noise_parameter_p <- paste0(stringr::str_sub(names(noise_parameters), 1, 1), noise_parameters, collapse = "")

    # create folder and file name
    folder_path <- glue::glue("00_sim_data/{dgp_setting$tau_variant}")
    file_name <- glue::glue('N{dgp_setting$N/1000}k.{dgp_setting$tau_variant}.X{stringr::str_sub(dgp_setting$x_dist, 1, 1)}.{dgp_setting$noise}_{noise_parameter_p}.rds')

    # check if file exists otherwise create DGP instructions ----
    if(file.exists(here::here(folder_path, file_name))){
        sim_settings <- readRDS(here::here(folder_path, file_name))
    }

    if(!file.exists(here::here(folder_path, file_name))){
        gen_df_sim_specs(
            n_dataset = dgp_setting$n_dataset, N = dgp_setting$N,
            D = dgp_setting$D,noise = dgp_setting$noise,
            noise_parameters = dgp_setting$noise_parameters[[1]],
            tau_variant = dgp_setting$tau_variant,
            sigma_variant = dgp_setting$sigma_variant,
            tau_effect = dgp_setting$tau_effect,
            x_dist = dgp_setting$x_dist)

        sim_settings <- readRDS(here::here(folder_path, file_name))
    }

    # simulate data
    data <- sim_settings %>%
        dplyr::mutate(sim_dataset = purrr::pmap(., simulate_data)) %>%
        # add information for algorithm
        dplyr::mutate(
            methods = list(methods),
            k = k,
            ntrees = ntrees)
    
    data %<>%
        dplyr::slice(1:5)

    # estimation
    data %<>%
        dplyr::mutate(result =
                          furrr::future_pmap(.,sim_est,
                                             seed = TRUE,
                                             .progress = TRUE,
                                             .options = furrr::furrr_options(
                                                 seed = TRUE,
                                                 conditions = "error",
                                                 packages = c("robustCausalTree", "rpart",
                                                              "partykit", "future", "mvtnorm",
                                                              "libcoin",
                                                              "dplyr", "tidyr", "magrittr")
                                             )))

    data %<>%
        dplyr::select(-sim_dataset)

    # save results
    save_folder <- glue::glue("02_results/01_sim/{dgp_setting$tau_variant}")
    # create folder if it does not exist
    dir.create(file.path(here::here(save_folder)), recursive = TRUE, showWarnings = FALSE)
    # save
    saveRDS(data, here::here(save_folder, file_name))
    
    if(verbose) {
        cli::cli_alert_success("Results successfully saved to \n {.path {here::here(save_folder, file_name)}}")
    }

}


sim_est <- function(sim_dataset, methods = methods, k = k, ntrees = ntrees, ...){

    on.exit(gc(), add = TRUE)  # clean on exit
    ## tries
    # sim_dataset <- data$sim_dataset[[1]]
    # methods <- data$methods[[1]]
    # k <- 4
    # ntress <- 4e3

    ###

    # get covariates
    X <- sim_dataset %>% dplyr::select(-c(Y, W, tau))

    # Sample ids for cross-fitting
    D = ncol(X) # number of covariats in X
    N = nrow(X) # number of obs. in dataset
    idx <- sample(c(1:k), N, replace = T) # index for k-fold cross validation

    # Centering ----

    center_dataset_Y <- cbind(Y=sim_dataset$Y, X)
    center_dataset_W <- cbind(W=sim_dataset$W, X)
    Y.hat <- rep(NaN, N)
    Y.res <- rep(NaN, N)
    W.hat <- rep(NaN, N)
    W.res <- rep(NaN, N)

    for(i in 1:k){ # centering for each fold

        Y.model <- grf::regression_forest(X = X[idx != i, ], Y = sim_dataset$Y[idx != i], num.trees = ntrees)
        Y.hat[idx == i] <- predict(Y.model, newdata = X[idx == i, ])$predictions
        Y.res[idx == i] <- center_dataset_Y$Y[idx == i] - Y.hat[idx == i]

        W.model <- grf::regression_forest(X = X[idx != i, ], Y = sim_dataset$W[idx != i], num.trees = ntrees)
        W.hat[idx == i] <- predict(W.model, newdata = X[idx == i, ])$predictions
        W.res[idx == i] <- center_dataset_W$W[idx == i] - W.hat[idx == i]
    }

    methods <- unlist(methods)

    # creating dataset
    df <- tibble::tibble(
        'sim_dataset' = list(sim_dataset),
        'Y.hat' = list(Y.hat),
        'W.hat' = list(W.hat),
        'Y.res' = list(Y.res),
        'W.res' = list(W.res),
        'k' = k,
        'idx' = list(idx),
        'ntrees' = ntrees,
        'split.R' = methods,
        'cv.opt' = methods,
        'ATE_estimation' = TRUE,
        'inference' = list(inference_options()) # for new function
    )


    df %<>%
        dplyr::mutate(purrr::pmap_df(.,est_and_eval, .progress = TRUE)) %>%
        dplyr::select(split.R, ends_with(c("cate", "ate")))

    return(df)
}


est_and_eval <- function(sim_dataset, Y.hat, W.hat, Y.res, W.res, k,
                         idx, ntrees, split.R, cv.opt, ATE_estimation,
                         inference, ...){

    on.exit(gc(), add = TRUE)  # clean on exit

    # tries
    # sim_dataset = df$sim_dataset[[1]]
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
    # inference = df$inference[[1]]
    # minsize_estim = 5

    # Estimation ----

    # Standard OLD version
    # estims <- robustCausalTree::robustCausalForest(Y=sim_dataset$Y,
    #                                                W=sim_dataset$W,
    #                                                Y.hat=Y.hat,
    #                                                W.hat=W.hat,
    #                                                Y.res=Y.res,
    #                                                W.res=W.res,
    #                                                X = sim_dataset %>% dplyr::select(-c(W, Y, tau)),
    #                                                k = k, idx = idx,
    #                                                ntrees = ntrees,
    #                                                split.R = split.R,
    #                                                cv.opt = cv.opt)

    # New function with BLB and IJ

    # robustCausalForest_BLB_IJ
    estims <- robustCausalForest_BLB_IJ_sim(Y=sim_dataset$Y,
                                        W=sim_dataset$W,
                                        Y.hat=Y.hat,
                                        W.hat=W.hat,
                                        Y.res=Y.res,
                                        W.res=W.res,
                                        X = sim_dataset %>%
                                            dplyr::select(-c(W, Y, tau)),
                                        k = k, idx = idx,
                                        ntrees = ntrees,
                                        split.R = split.R,
                                        cv.opt = cv.opt,
                                        .inference = inference)


    # Evaluation & Metric -----

    # coverage

    # conditional
    # names of standard errors
    cate_uncertainty <- list()

    std_err_names <- estims$pred_cate %>% select(starts_with("std.err_")) %>% names()

    for (i in seq_along(std_err_names)) {
        cate_uncertainty[[std_err_names[i]]] <- coverage_width_eval(
            estims$pred_cate$estimate,
            estims$pred_cate[[std_err_names[i]]],
            true_tau=sim_dataset$tau)
    }

    # average
    ate_uncertainty <- coverage_width_eval(estims$pred_ate["estimate"],
                                           estims$pred_ate["std.err"],
                                           true_tau = mean(sim_dataset$tau))

    # metrics
    df <- tibble::tibble(
        # saving split Rule
        'split.Rule' = split.R,
        # Conditional average treatment effect
        'cate' = list(estims$pred_cate$estimate),
        # Average conditional treatment effect
        'cate_avg' = mean(estims$pred_cate$estimate),
        # RMSE of conditional Average Treatment Effect
        "rmse_cate" = rmse_eval(estims$pred_cate$estimate, sim_dataset$tau),
        # abs. Bias of conditional Average Treatment Effect
        "absbias_cate" = absbias_eval(estims$pred_cate$estimate, sim_dataset$tau) %>% mean(),
        # CATE CI metrics
        "ci_cate" = list(cate_CI_metrics(cate_uncertainty)),
        # ATE ----
        "ate" = estims$pred_ate["estimate"],
        # Average prediction Error
        "err_ate" = estims$pred_ate["estimate"] - mean(sim_dataset$tau),
        # Coverage for ATE
        "coverage_ate" = ate_uncertainty$coverage,
        # CI width for ATE
        "width_ate" = ate_uncertainty$width
    )

    gc()  # clean

    return(df)

}


cate_CI_metrics <- function(cate_uncertainty, ...){
    cate_ci_metric <- tibble::enframe(cate_uncertainty) %>%
        tidyr::unnest(value) %>%
        dplyr::mutate(metric = names(value)) %>%
        tidyr::pivot_wider(names_from = metric, values_from = value) %>%
        group_by(name) %>%
        summarise(
            coverage_cate   = mean(unlist(coverage), na.rm = TRUE),
            avg_width_cate  = mean(unlist(width),  na.rm = TRUE),
            width_cate      = list(unlist(width))
        )

    return(cate_ci_metric)
}

