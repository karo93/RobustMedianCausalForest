source(here::here(
  "01_code", "02_simulation", "contamination",
  "archived_version", "01_code", "packages.R"
))

source(here::here(
  "01_code", "02_simulation", "contamination",
  "archived_version", "01_code", "helper.R"
))

package <- "Dev3"  # MSD

devtools::install(
  here::here(
    "01_code", "02_simulation", "contamination",
    "archived_version", "00_packages",
    paste0("causalTree", package)
  )
)

test_package <- causalTree::function_test()

data_dir <- here::here("data", "simulation", "contamination", "S1")
train_dir <- file.path(data_dir, "mean")
test_dir <- file.path(data_dir, "clean", "test")

result_dir <- here::here(
  "results", "simulation", "contamination", "S1", "mean", "MSD"
)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

dataset_names_test <- list.files(test_dir)[1]
dataset_names <- list.files(train_dir)
dataset_names <- dataset_names[grep("s2", dataset_names)]

set.seed(1)

splitRule <- "MSD"

for (k in 1:1) {
  for (j in 1:7) {
    dat1 <- readRDS(file.path(train_dir, dataset_names[j]))
    dat2 <- readRDS(file.path(test_dir, dataset_names_test))
    
    results <- lapply(1:max(dat1 %>% select(sim_id)), function(i) {
      print(i)
      
      datTrain <- dat1 %>%
        filter(sim_id == i) %>%
        select(2:ncol(dat1))
      
      datTest <- dat2 %>%
        filter(sim_id == i) %>%
        select(2:ncol(dat2))
      
      if (splitRule[k] == "CT") {
        preds <- CT_Dev(
          X = datTrain %>% select(X1:last_col()) %>% as.data.frame(),
          Y = datTrain$Yobs,
          T_ind = datTrain$T_ind,
          Xtest = datTest %>% select(X1:last_col()) %>% as.data.frame()
        )
      } else {
        preds <- user_Dev(
          X = datTrain %>% select(X1:last_col()) %>% as.data.frame(),
          Y = datTrain$Yobs,
          T_ind = datTrain$T_ind,
          Xtest = datTest %>% select(X1:last_col()) %>% as.data.frame()
        )
      }
      
      results <- preds %>%
        as_tibble() %>%
        mutate(
          true_CATE = datTest$CATE_true,
          true_ITE = datTest$Y1
        )
      
      return(results)
    })
    
    saveRDS(results, file.path(result_dir, dataset_names[j]))
  }
}

