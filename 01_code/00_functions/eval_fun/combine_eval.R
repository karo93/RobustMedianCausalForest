# function 
combine_evaluation <- function(analysis_path = "03_analysis/01_sim",
                               verbose = TRUE, ...){
  
  ## header, collect files
  
  absolut_path <- here::here(analysis_path)

  # check if folder exsists
  if (!dir.exists(absolut_path)) {
    cli::cli_abort(c(
      "Folder not found: {.path {absolut_path}}",
      "i" = "Please check that {.var analysis_path} is correct."
    ))
  }
  
  data <- tibble::tibble(
    rel_path = list.files(
      here::here(analysis_path),
      full.names = FALSE, recursive = TRUE,
      pattern = "*.Rdata")
    ) 
  
  # check if at least one case was already evaluated
  if (NROW(data) == 0) {
    cli::cli_abort(c(
      "No files found in {.path {here::here(analysis_path)}}",
      "x" = "No files matching pattern {.val *.Rdata} were found.",
      "i" = "Did you forget to run the analysis for this case?"
    ))
  }
  
  # parse information
  data %<>%
    # paste effect and noise from path
    dplyr::mutate(path_parts = stringr::str_split(rel_path, "/")) %>%
    tidyr::unnest_wider(path_parts, names_sep = "_") %>%
    dplyr::rename("effect" = path_parts_1,
                  "noise" = path_parts_2) %>%
    dplyr::select(-path_parts_3) 
  
  # load data
  data %<>%
    dplyr::mutate(
      data = purrr::map(rel_path, ~{
        e <- new.env()
        load(here::here(analysis_path, .x), envir = e)
        as.list(e)
        })
      )
  
  data %<>%
    ## unnest and drop path
    tidyr::unnest_wider(data) %>%
    tidyr::unnest_wider(metric) %>%
    dplyr::select(-rel_path) %>%
    # clean effect names to Scenarios
    dplyr::mutate(
      effect = case_when(
        effect == "complex_non_linear" ~ "S2",
        effect == "complex_linear" ~ "S1",
        effect == "no_effect" ~ "S3",
        effect == "cnl_interaction" ~ "S4",
        .default = effect),
      noise = case_when(
        noise == "heavy" ~ "Student-$t$",
        noise == "normal" ~ "Normal",
        noise == "laplace" ~ "Laplace",
        noise == "cauchy" ~ "Cauchy",
        .default = noise)
      ) %>%
    dplyr::mutate(effect = as.factor(effect))
  
  # CATE Analysis ---
  
  cate_analysis <-  combined_CATE(data, verbose)
  
  return(cate_analysis)

}


### CATE function 
combined_CATE <- function(data, verbose = TRUE, ...){
  
  # get data
  cate_eval <- data %>%
    dplyr::select(effect, noise, cate_eval) %>%
    tidyr::unnest(cate_eval) %>%
    tidyr::unnest(cate_eval) %>%
    tidyr::unnest(cate_eval)
  
  cate_table <- cate_eval %>%
    dplyr::select(effect, noise, split.R, starts_with("rmse_"), starts_with("absbias_")) %>%
    # combine Mean and SD for LaTeX Table
    tidyr::pivot_longer(
      cols = matches("_(mean|sd)$"),
      names_to = c("metric", "stat"),
      names_pattern = "^(.*)_(mean|sd)$",
      values_to = "val"
    ) %>%
    dplyr::mutate(val =
                    case_when(
                      stat == "sd" ~ sprintf("(%.2f)", val),
                      stat == "mean" ~ sprintf("%.2f", val)
                    )
    ) %>%
    tidyr::pivot_wider(names_from = metric, values_from = val) %>%
    dplyr::rename("MAE" = absbias,
                  "RMSE" = rmse) %>%
    tidyr::pivot_wider(
      id_cols = c(effect, noise, stat),
      names_from = split.R,
      values_from = c(RMSE, MAE),
      names_glue = "{split.R}.{.value}") %>%
    # sorting
    dplyr::select(effect, noise, stat,
                  starts_with("CT"),
                  starts_with("LMS"),
                  starts_with("MAD"),
                  starts_with("MSD")
                  ) %>%
    dplyr::arrange(effect, noise)
  
  ## Confidence intervals
  
  ci_eval <- cate_eval %>%
    dplyr::select(effect, noise, split.R, starts_with("BLB_"), starts_with("IJ_")) %>%
    tidyr::pivot_longer(
      cols = matches("_(mean|sd)$"),
      names_to = c("metric", "stat"),
      names_pattern = "^(.*)_(mean|sd)$",
      values_to = "val"
    ) %>%
    dplyr::mutate(val =
                    case_when(
                      stat == "sd" ~ sprintf("(%.2f)", val),
                      stat == "mean" ~ sprintf("%.2f", val))
    )  %>%
    tidyr::pivot_wider(names_from = metric, values_from = val) %>%
    tidyr::pivot_wider(
      id_cols = c(effect, noise, stat),
      names_from = split.R,
      values_from = c(
        BLB_uncor_coverage, BLB_uncor_ci_width,
        BLB_coverage, BLB_ci_width,
        IJ_uncor_coverage, IJ_uncor_ci_width,
        IJ_coverage, IJ_ci_width),
      names_glue = "{split.R}.{.value}") %>%
    # sorting
    dplyr::select(effect, noise, stat,
                  starts_with("CT"),
                  starts_with("LMS"),
                  starts_with("MAD"),
                  starts_with("MSD")) %>%
    dplyr::arrange(effect, noise)
    
  # just CT and LMS
  ci_eval_table_CT_LMS <- ci_eval %>%
    # sorting & selecting
    dplyr::select(effect, noise, stat,
                  starts_with("CT"),
                  starts_with("LMS")) %>%
    dplyr::arrange(effect, noise)
  
  ## MAD & MSD
  ci_eval_MAD_MSD <- ci_eval %>%
    # sorting & selectin
    dplyr::select(effect, noise, stat,
                  starts_with("MAD"),
                  starts_with("MSD")) %>%
    dplyr::arrange(effect, noise)
  
  return(
    list(
      "CATE Table"            = cate_table,
      "CI Table (CT & LMS)"   = ci_eval_table_CT_LMS,
      "CI Table (MAD & MSD)"  = ci_eval_MAD_MSD
      )
  )
  
}
