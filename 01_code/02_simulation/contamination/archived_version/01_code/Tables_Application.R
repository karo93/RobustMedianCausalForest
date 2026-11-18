library(dplyr)
library(purrr)
library(tibble)

archive_dir <- here::here(
  "01_code", "02_simulation", "contamination", "archived_version"
)

base_dir <- file.path(archive_dir, "02_raw_results", "Application")
output_dir <- file.path(
  archive_dir, "03_results_statistics", "Application", "ACIC"
)
#Change ACIC for other Application

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summarise_one_file <- function(file_path, dataset, contamination, method, p_c) {
  res_list <- readRDS(file_path)
  
  purrr::map_dfr(seq_along(res_list), function(i) {
    dat <- res_list[[i]]
    
    tibble(
      split = i,
      dataset = dataset,
      contamination = contamination,
      method = method,
      p_c = p_c,
      rmse_delta = dplyr::first(dat$rmse),
      mae_delta = dplyr::first(dat$mae),
      mape_delta = dplyr::first(dat$mape),
      correlation = dplyr::first(dat$correlation),
      mean_tau_contaminated = mean(dat$tau_hat, na.rm = TRUE),
      mean_tau_clean = mean(dat$tau_hat_clean, na.rm = TRUE),
      abs_delta_mean_tau = abs(mean_tau_contaminated - mean_tau_clean),
      sd_tau_contaminated = sd(dat$tau_hat, na.rm = TRUE),
      sd_tau_clean = sd(dat$tau_hat_clean, na.rm = TRUE),
      delta_sd_tau = sd_tau_contaminated - sd_tau_clean
    )
  })
}

make_application_contamination_table <- function(
    dataset,
    contamination,
    methods = c("CT", "LMS", "MAD", "MSD"),
    levels = c(5, 10, 15, 20, 25, 30, 35)) {
  
  out <- list()
  counter <- 1
  
  for (method in methods) {
    for (level in levels) {
      p_c <- level / 100
      
      file_path <- file.path(
        base_dir,
        dataset,
        contamination,
        method,
        paste0("results_variance_contaminated_", level, ".rds")
      )
      
      if (!file.exists(file_path)) {
        warning("Missing file: ", file_path)
        next
      }
      
      out[[counter]] <- summarise_one_file(
        file_path = file_path,
        dataset = dataset,
        contamination = contamination,
        method = method,
        p_c = p_c
      )
      
      counter <- counter + 1
    }
  }
  
  split_level <- dplyr::bind_rows(out)
  
  split_level %>%
    dplyr::group_by(dataset, contamination, method, p_c) %>%
    dplyr::summarise(
      rmse_delta = mean(rmse_delta, na.rm = TRUE),
      mae_delta = mean(mae_delta, na.rm = TRUE),
      mape_delta = mean(mape_delta, na.rm = TRUE),
      correlation = mean(correlation, na.rm = TRUE),
      abs_delta_mean_tau = mean(abs_delta_mean_tau, na.rm = TRUE),
      delta_sd_tau = mean(delta_sd_tau, na.rm = TRUE),
      n_splits = dplyr::n(),
      .groups = "drop"
    )
}

summary_acic_variance <- make_application_contamination_table(
  dataset = "ACIC",
  contamination = "variance_contamination"
)

saveRDS(
  summary_acic_variance,
  file.path(output_dir, "summary_acic_variance.rds")
)