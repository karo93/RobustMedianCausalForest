library(dplyr)
library(readr)
library(future)
library(future.apply)
library(rpart)

archive_dir <- here::here(
  "01_code", "02_simulation", "contamination", "archived_version"
)

source(file.path(archive_dir, "01_code", "packages.R"))
source(file.path(archive_dir, "01_code", "helper.R"))

package <- "Dev4"
method <- "LMS"

devtools::install(
  file.path(archive_dir, "00_packages", paste0("causalTree", package))
)

test_package <- causalTree::function_test()

data_dir <- here::here("00_data", "00_Application", "ACIC")
clean_dir <- file.path(data_dir, "clean_data")
contam_dir <- file.path(data_dir, "variance_contamination")

results_dir <- file.path(
  archive_dir, "02_raw_results", "Application", "ACIC"
)

clean_results_dir <- file.path(results_dir, "clean_results", method)
contam_results_dir <- file.path(results_dir, "variance_contamination", method)

dir.create(clean_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(contam_results_dir, recursive = TRUE, showWarnings = FALSE)

splittingrule <- "user"
p_values <- c(5, 10, 15, 20, 25, 30, 35)
n_splits <- 100

ACIC_clean <- read_csv(file.path(clean_dir, "ACIC_clean.csv"))
colnames(ACIC_clean) <- make.names(colnames(ACIC_clean))

feature_cols_prop <- setdiff(colnames(ACIC_clean), c("Y", "W"))
propensity_formula <- as.formula(
  paste("W ~", paste(feature_cols_prop, collapse = " + "))
)

propensity_model <- glm(
  propensity_formula,
  data = ACIC_clean,
  family = binomial
)

clean_partition <- readRDS(file.path(clean_dir, "clean.rds"))

future::plan(multisession, workers = 9)

clean_tau_hat_list <- future_lapply(seq_len(n_splits), function(i) {
  library(rpart)
  
  cat("Training clean model on split", i, "\n")
  
  datTrain <- clean_partition[[i]]$train
  datTest <- clean_partition[[i]]$test
  
  colnames(datTrain) <- make.names(colnames(datTrain))
  colnames(datTest) <- make.names(colnames(datTest))
  
  feature_cols <- setdiff(colnames(datTrain), c("Y", "W"))
  f_2 <- paste(feature_cols, collapse = " + ")
  
  Ttrain <- datTrain$W
  
  forest_clean <- causalTree::causalForest(
    as.formula(paste("Y ~", f_2)),
    data = datTrain,
    treatment = Ttrain,
    split.Rule = splittingrule,
    split.Honest = TRUE,
    split.Bucket = FALSE,
    bucketNum = 5,
    bucketMax = 100,
    cv.option = splittingrule,
    cv.Honest = TRUE,
    minsize = 2L,
    split.alpha = 0.5,
    cv.alpha = 0.5,
    sample.size.total = floor(nrow(datTrain) / 2),
    sample.size.train.frac = 0.5,
    mtry = ceiling(ncol(datTrain) / 3),
    nodesize = 10,
    num.trees = 1000,
    ncolx = length(feature_cols),
    ncov_sample = floor(length(feature_cols) / 3)
  )
  
  Xtest <- datTest[, feature_cols]
  as.numeric(predict(forest_clean, Xtest))
}, future.seed = TRUE)

saveRDS(
  clean_tau_hat_list,
  file.path(clean_results_dir, "tau_hat_clean_list.rds")
)

tau_hat_clean_list <- readRDS(
  file.path(clean_results_dir, "tau_hat_clean_list.rds")
)

for (p in p_values) {
  contaminated_splits <- readRDS(
    file.path(contam_dir, paste0("variance_contaminated_", p, ".rds"))
  )
  
  cat("Starting contamination level p =", p, "\n")
  
  results_all <- future_lapply(seq_len(n_splits), function(i) {
    library(rpart)
    
    datTrain <- contaminated_splits[[i]]$train
    datTest <- contaminated_splits[[i]]$test
    
    colnames(datTrain) <- make.names(colnames(datTrain))
    colnames(datTest) <- make.names(colnames(datTest))
    
    feature_cols <- setdiff(colnames(datTrain), c("Y", "W"))
    f_2 <- paste(feature_cols, collapse = " + ")
    
    Ttrain <- datTrain$W
    Xtrain <- datTrain[, feature_cols]
    Xtest <- datTest[, feature_cols]
    
    forest_honest <- causalTree::causalForest(
      as.formula(paste("Y ~", f_2)),
      data = datTrain,
      treatment = Ttrain,
      split.Rule = splittingrule,
      split.Honest = TRUE,
      split.Bucket = FALSE,
      bucketNum = 5,
      bucketMax = 100,
      cv.option = splittingrule,
      cv.Honest = TRUE,
      minsize = 2L,
      split.alpha = 0.5,
      cv.alpha = 0.5,
      sample.size.total = floor(nrow(datTrain) / 2),
      sample.size.train.frac = 0.5,
      mtry = ceiling(ncol(datTrain) / 3),
      nodesize = 10,
      num.trees = 1000,
      ncolx = length(Xtrain),
      ncov_sample = floor(length(Xtrain) / 3)
    )
    
    tau_hat <- as.numeric(predict(forest_honest, Xtest))
    
    propens <- predict(
      propensity_model,
      newdata = datTest,
      type = "response"
    )
    propens <- pmin(pmax(propens, 1e-6), 1 - 1e-6)
    
    ystar <- datTest$Y *
      (datTest$W - propens) / (propens * (1 - propens))
    
    tau_hat_clean <- tau_hat_clean_list[[i]]
    
    mape <- mean(
      abs((tau_hat - tau_hat_clean) / (abs(tau_hat_clean) + 1e-6))
    ) * 100
    
    rmse <- sqrt(mean((tau_hat - tau_hat_clean)^2, na.rm = TRUE))
    mae <- mean(abs(tau_hat - tau_hat_clean))
    cor_tau <- cor(tau_hat, tau_hat_clean)
    
    tibble(
      tau_hat = tau_hat,
      tau_hat_clean = tau_hat_clean,
      ystar = ystar,
      mape = mape,
      rmse = rmse,
      mae = mae,
      correlation = cor_tau
    )
  }, future.seed = TRUE)
  
  saveRDS(
    results_all,
    file.path(
      contam_results_dir,
      paste0("results_variance_contaminated_", p, ".rds")
    )
  )
}

future::plan(sequential)