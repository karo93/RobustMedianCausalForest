# wrapper function
analysis_fun <- function(dgp_object, verbose = TRUE, ...){
    
    # Parameters
    tau_variant <- dgp_object$tau_variant
    noise <- dgp_object$noise
    noise_parameters <- dgp_object$noise_parameters[[1]]
    N <- dgp_object$N 
    D <- dgp_object$D
    tau_effect <- dgp_object$tau_effect
    sigma_variant <- dgp_object$sigma_variant
    x_dist <- dgp_object$x_dist
    
    # prepare noise parameter for naming
    noise_parameter_p <- paste0(stringr::str_sub(names(noise_parameters), 1, 1), noise_parameters, collapse = "")

    # create folder and file name
    folder_path <- glue::glue("02_results/01_sim/{tau_variant}")
    file_name <- glue::glue('N{N/1000}k.{tau_variant}.X{stringr::str_sub(x_dist, 1, 1)}.{noise}_{noise_parameter_p}')
    save_path <- glue::glue("03_analysis/01_sim/{tau_variant}/{noise}")
    path_file_name <- here::here(glue::glue("{save_path}/{file_name}.Rdata"))
    
    dir.create(file.path(save_path), recursive = TRUE, showWarnings = FALSE)

    # load RDS sim data
    data_set <- load_rds_safely(folder_path, file_name)
#    data_set <- readRDS(here::here(folder_path, glue::glue("{file_name}.Rds")))

    # data wrangling ----
    results <- data_set %>%
        dplyr::select(dataset_num, tau_variant, result, noise) %>%
        tidyr::unnest(result)

    # get ci_data
    ci_results <- results %>%
        tidyr::unnest(ci_cate) %>%
        # clean method names
        dplyr::rename('CI_method' = name) %>%
        dplyr::mutate(CI_method = stringr::str_remove(CI_method, "std.err_"))

    # average effects analysis ----
    ate_analysis <- ate_analysis(results)

    # conditional effects analysis ----
    cate_analysis <- cate_analysis(ci_results)

    metric <- tibble::tibble(
        "ate_eval" = list(ate_analysis),
        "cate_eval" = list(cate_analysis)
        )
    
    save(metric, file = path_file_name)
    
    if (verbose) {
        cli::cli_alert_success("Evaluation successfully saved to \n {.path {path_file_name}}")
    }
    
    return(metric)
}

## ate function
ate_analysis <- function(results, ...){

    # score function summary ----
    ate_summary <- results %>%
        dplyr::select(dataset_num, split.R, ate, err_ate, coverage_ate, width_ate) %>%
        dplyr::distinct(dataset_num, split.R, .keep_all = TRUE) %>%
        dplyr::rename(error = err_ate,
                      coverage = coverage_ate,
                      ci_width = width_ate) %>%
        dplyr::summarise(
            across(
                c(ate, error, coverage, ci_width),
                list(
                    mean   = function(x) mean(x, na.rm = TRUE),
                    sd     = function(x) sd(x, na.rm = TRUE)
                )
            ),
            .by = c(split.R)
        )
    
    # return objekt
    #return(tibble::tibble("ate_eval" = list(ate_summary)))
    return(ate_summary)
}


## CATE Analysis
cate_analysis <- function(ci_results, save_path, file_name, render_plot = FALSE,
                          plot_parameters = plot_parameters, ...){

    # score function summary
    cate_summary <- ci_results %>%
        dplyr::select(dataset_num, split.R, rmse_cate, absbias_cate) %>%
        # delete duplicates due too two CI methods (do not playe a role here)
        dplyr::distinct(dataset_num, split.R, .keep_all = TRUE) %>%
        dplyr::rename(rmse = rmse_cate,
                      absbias = absbias_cate) %>%
        dplyr::summarise(
            across(
                c(rmse, absbias),
                list(
                    mean   = function(x) mean(x, na.rm = TRUE),
                    sd     = function(x) sd(x, na.rm = TRUE)
                )
            ),
            .by = c(split.R)
        )

    # summary stats CI
    ci_summary <- ci_results %>%
        dplyr::select(split.R, CI_method, coverage_cate, avg_width_cate) %>%
        dplyr::rename(coverage = coverage_cate,
                      ci_width = avg_width_cate) %>%
        dplyr::summarise(
            across(
                c(coverage, ci_width),
                list(
                    mean   = function(x) mean(x, na.rm = TRUE),
                    sd     = function(x) sd(x, na.rm = TRUE)
                )
            ),
            .by = c(split.R, CI_method)
        ) %>%
        tidyr::pivot_wider(names_from = CI_method, values_from = coverage_mean:ci_width_sd,  names_glue = "{CI_method}_{.value}") %>%
        select(
            split.R,
            starts_with("BLB_"),
            starts_with("IJ_")
        )

    # combine score results and CI results
    cate_summary %<>%
        dplyr::left_join(ci_summary,
                         by = "split.R")

   
    return(tibble::tibble("cate_eval" = list(cate_summary)))

}
