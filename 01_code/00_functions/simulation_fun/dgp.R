gen_df_sim_specs <- function(N, D, noise, noise_parameters,
                             tau_variant, sigma_variant = 2, # sigmund variant
                             tau_effect = "het", x_dist =  'normal', n_dataset = 100) {

    # prepare noise parameter for naming
    noise_parameter_p <- paste0(stringr::str_sub(names(noise_parameters), 1, 1), noise_parameters, collapse = "")
    # create folder if not present
    folder_path <- here::here('00_sim_data', tau_variant)

    dir.create(file.path(folder_path), recursive = TRUE, showWarnings = FALSE)

    df <- tibble::tibble(
        n = N,
        d = D,
        noise = noise,
        noise_parameters = list(noise_parameters),
        tau_effect = tau_effect,
        tau_variant = tau_variant,
        sigma_variant = 2, # sigmund variant
        x_dist =  x_dist,
        dataset_num = list(1:n_dataset)) %>%
        tidyr::unnest(dataset_num) %>%
        dplyr::mutate(seed = sample(1e3:1e4-1, dplyr::n(), replace = TRUE))

    # save df
    saveRDS(df, glue::glue('{folder_path}/N{N/1000}k.{tau_variant}.X{stringr::str_sub(x_dist, 1, 1)}.{noise}_{noise_parameter_p}.rds'))
}
