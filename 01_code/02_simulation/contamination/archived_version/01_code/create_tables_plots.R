library(knitr)
library(kableExtra)

methods_vec <- c("CT", "LMS", "MAD", "MSD")
contam_levels_vec <- c(0, 5, 10, 15, 20, 25, 30, 35)

base_sim_path <- here::here("02_raw_results")
output_dir <- here::here("03_results_statistics")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scenario_files <- c(
  S1 = "sim_low_s2_train",
  S2 = "sim_low_s3_train",
  S3 = "sim_low_s4_train",
  S4 = "sim_s3_train"
)

read_simulation_results <- function(
    scenario, contamination_folder, contamination_label,
    file_stub, methods_use = methods_vec,
    contam_levels_use = contam_levels_vec
) {
  
  map_dfr(methods_use, function(method_i) {
    map_dfr(contam_levels_use, function(p_i) {
      
      folder <- if (p_i == 0) "00_Clean" else contamination_folder
      filename <- paste0(file_stub, if (p_i == 0) "" else paste0("_", p_i), ".rds")
      
      file_path <- file.path(
        base_sim_path, scenario, folder, method_i, filename
      )
      
      if (!file.exists(file_path)) {
        warning("Missing file: ", file_path)
        return(tibble())
      }
      
      results <- purrr::compact(readRDS(file_path))
      
      map_dfr(seq_along(results), function(rep_id) {
        df <- results[[rep_id]]
        
        if (!all(c("true_CATE", "value") %in% names(df))) {
          stop("Missing prediction or true CATE column in: ", file_path)
        }
        
        df %>%
          summarise(
            rmse = sqrt(mean((value - true_CATE)^2, na.rm = TRUE)),
            mae = mean(abs(value - true_CATE), na.rm = TRUE),
            mape = mean(
              abs((value - true_CATE) / pmax(abs(true_CATE), 1)),
              na.rm = TRUE
            ) * 100
          ) %>%
          mutate(
            scenario = scenario,
            contamination = contamination_label,
            p_c = p_i,
            method = method_i,
            rep = rep_id
          )
      })
    })
  })
}

sim_results_rep <- map_dfr(names(scenario_files), function(scenario) {
  bind_rows(
    read_simulation_results(
      scenario, "01_Mean", "Mean-shift", scenario_files[[scenario]]
    ),
    read_simulation_results(
      scenario, "02_Variance", "Variance", scenario_files[[scenario]]
    )
  )
})

sim_results_summary <- sim_results_rep %>%
  group_by(scenario, contamination, method, p_c) %>%
  summarise(
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    sd_mae = sd(mae, na.rm = TRUE),
    mean_mape = mean(mape, na.rm = TRUE),
    sd_mape = sd(mape, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(method, levels = methods_vec),
    contamination = factor(contamination, levels = c("Mean-shift", "Variance")),
    scenario = factor(scenario, levels = names(scenario_files))
  )

sim_results_pct <- sim_results_summary %>%
  group_by(scenario, contamination, method) %>%
  mutate(
    rmse_clean = mean_rmse[p_c == 0][1],
    mae_clean = mean_mae[p_c == 0][1],
    mape_clean = mean_mape[p_c == 0][1],
    prmse = 100 * (mean_rmse - rmse_clean) / rmse_clean,
    pmae = 100 * (mean_mae - mae_clean) / mae_clean,
    pmape = 100 * (mean_mape - mape_clean) / mape_clean
  ) %>%
  ungroup()

method_cols <- c(
  CT = "#b2182b",
  LMS = "#ef8a62",
  MAD = "#91bfdb",
  MSD = "#2166ac"
)

sim_raw_long <- sim_results_summary %>%
  select(scenario, contamination, method, p_c,
         mean_mae, mean_mape, mean_rmse) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "metric",
    names_prefix = "mean_",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric, mae = "MAE", mape = "MAPE", rmse = "RMSE"),
    metric = factor(metric, levels = c("MAE", "MAPE", "RMSE"))
  )

sim_pct_long <- sim_results_pct %>%
  select(scenario, contamination, method, p_c, pmae, pmape, prmse) %>%
  pivot_longer(
    cols = c(pmae, pmape, prmse),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric, pmae = "PMAE", pmape = "PMAPE", prmse = "PRMSE"
    ),
    metric = factor(metric, levels = c("PMAE", "PMAPE", "PRMSE"))
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

make_sim_raw_plot <- function(scenario_name) {
  p <- sim_raw_long %>%
    filter(scenario == scenario_name) %>%
    ggplot(aes(x = p_c, y = value, color = method, group = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_manual(values = method_cols) +
    facet_grid(metric ~ contamination, scales = "free_y") +
    labs(
      x = "Contamination Level (%)",
      y = "Metric Value",
      color = "Method"
    ) +
    plot_theme
  
  ggsave(
    file.path(output_dir, paste0("Simulation_Raw_", scenario_name, ".png")),
    p, width = 10, height = 6, dpi = 300
  )
  
  p
}

make_sim_pct_plot <- function(scenario_name) {
  p <- sim_pct_long %>%
    filter(scenario == scenario_name) %>%
    ggplot(aes(x = p_c, y = value, color = method, group = method)) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_manual(values = method_cols) +
    facet_grid(metric ~ contamination, scales = "free_y") +
    labs(
      x = "Contamination Level (%)",
      y = "Increase (%) relative to clean baseline",
      color = "Method"
    ) +
    plot_theme
  
  ggsave(
    file.path(output_dir, paste0("Simulation_Percentage_", scenario_name, ".png")),
    p, width = 10, height = 6, dpi = 300
  )
  
  p
}

for (scenario in names(scenario_files)) {
  make_sim_raw_plot(scenario)
  make_sim_pct_plot(scenario)
}

sim_clean_table <- sim_results_summary %>%
  filter(p_c == 0) %>%
  select(scenario, method, mean_rmse, mean_mae, mean_mape) %>%
  mutate(across(c(mean_rmse, mean_mae, mean_mape), ~ round(.x, 3))) %>%
  arrange(scenario, method)

sim_clean_latex <- sim_clean_table %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    col.names = c(
      "Scenario", "Method",
      "$\\overline{\\operatorname{RMSE}}$",
      "$\\overline{\\operatorname{MAE}}$",
      "$\\overline{\\operatorname{MAPE}}$"
    ),
    caption = "Clean baseline performance in the simulation study.",
    label = "tab:sim_clean_baseline"
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 9
  )

writeLines(sim_clean_latex, file.path(output_dir, "table_sim_clean_baseline.tex"))

sim_raw_table <- sim_results_summary %>%
  select(scenario, contamination, method, p_c,
         mean_rmse, mean_mae, mean_mape) %>%
  mutate(across(c(mean_rmse, mean_mae, mean_mape), ~ round(.x, 3))) %>%
  arrange(scenario, contamination, method, p_c)

sim_raw_latex <- sim_raw_table %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    col.names = c(
      "Scenario", "Contamination", "Method", "$p_c$",
      "$\\overline{\\operatorname{RMSE}}$",
      "$\\overline{\\operatorname{MAE}}$",
      "$\\overline{\\operatorname{MAPE}}$"
    ),
    caption = "Raw simulation performance under selected contamination levels.",
    label = "tab:sim_raw_selected"
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 8
  )

writeLines(sim_raw_latex, file.path(output_dir, "table_sim_raw.tex"))

sim_pct_table <- sim_results_pct %>%
  select(scenario, contamination, method, p_c, prmse, pmae, pmape) %>%
  mutate(across(c(prmse, pmae, pmape), ~ round(.x, 3))) %>%
  arrange(scenario, contamination, method, p_c)

sim_pct_latex <- sim_pct_table %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    col.names = c(
      "Scenario", "Contamination", "Method", "$p_c$",
      "$\\operatorname{PRMSE}$",
      "$\\operatorname{PMAE}$",
      "$\\operatorname{PMAPE}$"
    ),
    caption = "Percentage increase in simulation error metrics relative to the clean baseline.",
    label = "tab:sim_pct_selected"
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 8
  )

writeLines(sim_pct_latex, file.path(output_dir, "table_sim_pct.tex"))