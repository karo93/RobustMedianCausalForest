methods <- c("CT", "LMS", "MAD", "MSD")
contam_levels <- c(5, 10, 15, 20, 25, 30, 35)
datasets <- c("ACIC", "Lalonde", "BreastCancer")

contamination_specs <- tibble(
  contamination_folder = c("mean_contamination", "variance_contamination"),
  type = c("Mean-shift", "Variance"),
  file_prefix = c("results_mean_contaminated", "results_variance_contaminated")
)

archive_dir <- here::here(
  "01_code", "02_simulation", "contamination", "archived_version"
)

base_path <- file.path(archive_dir, "02_raw_results", "Application")
output_dir <- file.path(archive_dir, "03_results_statistics", "Application")

method_cols <- c(
  "CT"  = "#b2182b",
  "LMS" = "#ef8a62",
  "MAD" = "#91bfdb",
  "MSD" = "#2166ac"
)

plot_theme <- theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "#004c93", color = "#004c93"),
    strip.text = element_text(face = "bold", color = "white"),
    panel.border = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey92"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.key = element_blank()
  )
read_application_results <- function(dataset_name,
                                     contamination_folder,
                                     type_label,
                                     file_prefix,
                                     methods_use = methods,
                                     contam_levels_use = contam_levels,
                                     base_path_use = base_path) {
  
  map_dfr(methods_use, function(method_i) {
    map_dfr(contam_levels_use, function(p_i) {
      
      file_path <- file.path(
        base_path_use,
        dataset_name,
        contamination_folder,
        method_i,
        paste0(file_prefix, "_", p_i, ".rds")
      )
      
      if (!file.exists(file_path)) {
        warning("Missing file: ", file_path)
        return(tibble())
      }
      
      results <- readRDS(file_path)
      results <- compact(results)
      
      map_dfr(seq_along(results), function(rep_id) {
        df <- results[[rep_id]]
        
        if (is.null(df)) return(tibble())
        
        df %>%
          summarise(
            mape = mean(abs((tau_hat - tau_hat_clean) / pmax(abs(tau_hat_clean), 1)),
                        na.rm = TRUE) * 100,
            rmse = sqrt(mean((tau_hat - tau_hat_clean)^2, na.rm = TRUE)),
            mae  = mean(abs(tau_hat - tau_hat_clean), na.rm = TRUE),
            correlation = cor(tau_hat, tau_hat_clean, use = "complete.obs")
          ) %>%
          mutate(
            dataset = dataset_name,
            type = type_label,
            p = p_i,
            method = method_i,
            rep = rep_id
          )
      })
    })
  })
}
all_results_rep <- map_dfr(datasets, function(dataset_i) {
  pmap_dfr(
    contamination_specs,
    function(contamination_folder, type, file_prefix) {
      read_application_results(
        dataset_name = dataset_i,
        contamination_folder = contamination_folder,
        type_label = type,
        file_prefix = file_prefix
      )
    }
  )
})
all_results_summary <- all_results_rep %>%
  group_by(dataset, type, method, p) %>%
  summarise(
    mean_mape = mean(mape, na.rm = TRUE),
    sd_mape   = sd(mape, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse   = sd(rmse, na.rm = TRUE),
    mean_mae  = mean(mae, na.rm = TRUE),
    sd_mae    = sd(mae, na.rm = TRUE),
    mean_corr = mean(correlation, na.rm = TRUE),
    sd_corr   = sd(correlation, na.rm = TRUE),
    n_reps    = n(),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(method, levels = c("CT", "LMS", "MAD", "MSD")),
    type = factor(type, levels = c("Mean-shift", "Variance")),
    dataset = factor(dataset, levels = c("ACIC", "Lalonde", "BreastCancer"))
  )
all_results_long <- all_results_summary %>%
  select(dataset, type, method, p, mean_mape, mean_rmse, mean_mae) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "metric",
    names_prefix = "mean_",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      mae  = "MAE_delta",
      mape = "MAPE_delta",
      rmse = "RMSE_delta"
    ),
    metric = factor(metric, levels = c("MAE_delta", "MAPE_delta", "RMSE_delta")),
    method = factor(method, levels = c("CT", "LMS", "MAD", "MSD")),
    type = factor(type, levels = c("Mean-shift", "Variance"))
  )
make_application_plot <- function(dataset_name) {
  
  plot_data <- all_results_long %>%
    filter(dataset == dataset_name)
  
  p <- ggplot(plot_data, aes(x = p, y = value, color = method, group = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_manual(values = method_cols) +
    facet_grid(
      metric ~ type,
      scales = "free_y",
      labeller = labeller(
        metric = as_labeller(c(
          MAE_delta  = "MAE[Delta]",
          MAPE_delta = "MAPE[Delta]",
          RMSE_delta = "RMSE[Delta]"
        ), label_parsed)
      )
    ) +
    labs(
      x = "Contamination Level (%)",
      y = "Metric Value",
      color = "Method"
    ) +
    plot_theme
  
  dataset_dir <- file.path(output_dir, dataset_name)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    file.path(dataset_dir, paste0(dataset_name, "_Robust.png")),
    p, width = 10, height = 6, dpi = 300
  )
  
  p
}
plots <- datasets %>%
  set_names(datasets) %>%
  map(make_application_plot)

ACIC_Robust <- plots$ACIC
Lalonde_Robust <- plots$Lalonde
Breastcancer_Robust <- plots$BreastCancer

###################### plots relative to best
relative_to_best <- all_results_summary %>%
  select(dataset, type, method, p,
         mean_mae, mean_mape, mean_rmse) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(dataset, type, p, metric) %>%
  mutate(
    best_value = min(value, na.rm = TRUE),
    pct_increase_best = 100 * (value - best_value) / best_value
  ) %>%
  ungroup() %>%
  mutate(
    metric = recode(
      metric,
      mean_mae  = "MAE\u0394",
      mean_mape = "MAPE\u0394",
      mean_rmse = "RMSE\u0394"
    ),
    metric = factor(metric, levels = c("MAE\u0394", "MAPE\u0394", "RMSE\u0394")),
    method = factor(method, levels = c("CT", "LMS", "MAD", "MSD")),
    dataset = factor(dataset, levels = c("Lalonde", "BreastCancer", "ACIC")),
    type = factor(type, levels = c("Mean-shift", "Variance"))
  )

make_relative_boxplot_dataset <- function(dataset_name, output_file) {
  
  p <- relative_to_best %>%
    filter(dataset == dataset_name) %>%
    ggplot(aes(x = method, y = pct_increase_best, color = method)) +
    geom_boxplot(width = 0.6, outlier.alpha = 0.7) +
    scale_color_manual(values = method_cols) +
    facet_grid(metric ~ type, scales = "free_y") +
    labs(
      x = "Method",
      y = "Increase (%) relative to best",
      color = "Method"
    ) +
    theme_minimal(base_size = 13) +
    theme(    strip.background = element_rect(fill = "#004c93", color = "#004c93"),
              strip.text = element_text(face = "bold", color = "white"),
              legend.position = "none",
              panel.grid.minor = element_blank()
    )
  
  ggsave(output_file, p, width = 9, height = 6, dpi = 300)
  p
}
p_rel_lalonde <- make_relative_boxplot_dataset(
  "Lalonde",
  "Relative_Best_Lalonde.png"
)

p_rel_bc <- make_relative_boxplot_dataset(
  "BreastCancer",
  "Relative_Best_BreastCancer.png"
)

p_rel_acic <- make_relative_boxplot_dataset(
  "ACIC",
  "Relative_Best_ACIC.png"
)
