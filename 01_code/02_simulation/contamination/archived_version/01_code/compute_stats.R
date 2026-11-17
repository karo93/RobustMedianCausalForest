scenario <- "S4"
contamination_type <- "00_Clean"
splitrule_names <- c("CT", "LMS", "MAD", "MSD")

input_dir <- here::here("02_raw_results", scenario, contamination_type)
output_dir <- here::here("03_results_statistics", scenario, contamination_type)

dataset_names <- list.files(file.path(input_dir, "CT"))

for (i in seq_along(splitrule_names)) {
  splitrule <- splitrule_names[i]
  
  result_dir <- file.path(output_dir, splitrule)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  lapply(dataset_names, function(k) {
    dataset <- readRDS(file.path(input_dir, splitrule, k))
    
    dplyr::bind_rows(dataset) %>%
      mutate(id = rep(seq_along(dataset), each = nrow(dataset[[1]]))) %>%
      group_by(id) %>%
      summarise(
        precision_RMSE = caret::RMSE(value, true_CATE),
        precision_MAPE = mean(abs(value - true_CATE))
      ) %>%
      saveRDS(file = file.path(result_dir, k))
  })
}