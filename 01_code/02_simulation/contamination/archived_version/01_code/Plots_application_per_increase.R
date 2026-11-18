library(tidyverse)

archive_dir <- here::here(
  "01_code", "02_simulation", "contamination", "archived_version"
)

base_path <- file.path(archive_dir, "02_raw_results", "Application")
output_dir <- file.path(archive_dir, "03_results_statistics", "Application")

methods <- c("CT", "MAD", "MSD", "LMS")
contam_levels <- c(5, 10, 15, 20, 25, 30, 35)
datasets <- c("ACIC", "Lalonde", "BreastCancer")

contamination_specs <- tibble(
  contamination_folder = c("mean_contamination", "variance_contamination"),
  type = c("Mean", "Variance"),
  file_prefix = c(
    "results_mean_contaminated",
    "results_variance_contaminated"
  )
)

colors <- c(
  CT = "#b2182b",
  LMS = "#ef8a62",
  MAD = "#91bfdb",
  MSD = "#2166ac"
)

summarise_results <- function(dataset, contamination, type, file_prefix) {
  map_dfr(methods, function(method) {
    map_dfr(contam_levels, function(p) {
      file_path <- file.path(
        base_path,
        dataset,
        contamination,
        method,
        paste0(file_prefix, "_", p, ".rds")
      )
      
      if (!file.exists(file_path)) {
        warning("Missing file: ", file_path)
        return(tibble())
      }
      
      results <- purrr::compact(readRDS(file_path))
      
      map_dfr(seq_along(results), function(rep_id) {
        df <- results[[rep_id]]
        
        df %>%
          summarise(
            mape = mean(
              abs((tau_hat - tau_hat_clean) /
                    pmax(abs(tau_hat_clean), 1))
            ) * 100,
            rmse = sqrt(
              mean((tau_hat - tau_hat_clean)^2, na.rm = TRUE)
            ),
            mae = mean(
              abs(tau_hat - tau_hat_clean), na.rm = TRUE
            ),
            correlation = cor(
              tau_hat, tau_hat_clean, use = "complete.obs"
            )
          ) %>%
          mutate(
            p = p,
            method = method,
            rep = rep_id,
            dataset = dataset,
            type = type
          )
      })
    })
  })
}

all_results <- pmap_dfr(
  expand_grid(
    dataset = datasets,
    contamination_specs
  ),
  function(dataset, contamination_folder, type, file_prefix) {
    summarise_results(
      dataset = dataset,
      contamination = contamination_folder,
      type = type,
      file_prefix = file_prefix
    )
  }
)

summary_results <- all_results %>%
  group_by(dataset, type, method, p) %>%
  summarise(
    mean_mape = mean(mape, na.rm = TRUE),
    sd_mape = sd(mape, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_corr = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )

all_results_long <- summary_results %>%
  select(dataset, type, method, p, mean_mape, mean_rmse, mean_mae) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "metric",
    names_prefix = "mean_",
    values_to = "value"
  )

for (dataset in datasets) {
  plot_data <- all_results_long %>%
    filter(.data$dataset == .env$dataset)
  
  robustness_plot <- ggplot(
    plot_data,
    aes(x = p, y = value, color = method)
  ) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colors) +
    facet_grid(metric ~ type, scales = "free_y") +
    labs(
      x = "Contamination Level (%)",
      y = "Metric Value",
      color = "Method"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.background = element_rect(
        fill = "#004c93", color = "#004c93"
      ),
      strip.text = element_text(
        face = "bold", color = "white"
      ),
      panel.border = element_blank(),
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_line(color = "grey92"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.key = element_blank()
    )
  
  dataset_dir <- file.path(output_dir, dataset)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename = file.path(dataset_dir, paste0(dataset, "_Robust.png")),
    plot = robustness_plot,
    width = 10,
    height = 6,
    dpi = 300
  )
}