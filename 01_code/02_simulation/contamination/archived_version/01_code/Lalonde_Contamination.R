set.seed(123)

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

data_dir <- here::here("00_data", "00_Application", "Lalonde")
clean_dir <- file.path(data_dir, "clean_data")
contam_dir <- file.path(data_dir, "mean_contamination")

results_dir <- file.path(
  archive_dir, "02_raw_results", "Application", "Lalonde"
)

clean_results_dir <- file.path(results_dir, "clean_results", method)
contam_results_dir <- file.path(results_dir, "mean_contamination", method)

dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(contam_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clean_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(contam_results_dir, recursive = TRUE, showWarnings = FALSE)

lalonde_data <- lalonde::nsw[, -1]

n_splits <- 100
test_fraction <- 0.25
n <- nrow(lalonde_data)

splits_list <- vector("list", n_splits)

for (i in 1:n_splits) {
  test_id <- sample(1:n, size = floor(test_fraction * n))
  train_id <- setdiff(1:n, test_id)
  
  splits_list[[i]] <- list(
    train = as.data.frame(lalonde_data[train_id, ]),
    test = as.data.frame(lalonde_data[test_id, ])
  )
}

saveRDS(splits_list, file.path(clean_dir, "clean.rds"))

data_contamination <- function(data, column, p, newValues) {
  numReplace <- ceiling(nrow(data) * p / 100)
  replaceIndices <- sample(1:nrow(data), size = numReplace)
  
  data[replaceIndices, column] <- sample(
    newValues, size = numReplace, replace = TRUE
  )
  
  return(data)
}

p_values <- c(5, 10, 15, 20, 25, 30, 35)

for (p in p_values) {
  contaminated_splits <- vector("list", length(splits_list))
  
  for (i in seq_along(splits_list)) {
    train_data <- splits_list[[i]]$train
    test_data <- splits_list[[i]]$test
    
    sigma <- sd(train_data$re78)
    
    Ycontamin <- train_data$re78 +
      rnorm(nrow(train_data), mean = 5 * sigma, sd = sigma)
    
    # Variance contamination
    # Ycontamin <- train_data$re78 +
    #   rnorm(nrow(train_data), mean = 0, sd = 5 * sigma)
    
    train_contaminated <- data_contamination(
      train_data, column = "re78", p = p, newValues = Ycontamin
    )
    
    contaminated_splits[[i]] <- list(
      train = train_contaminated,
      test = test_data
    )
  }
  
  saveRDS(
    contaminated_splits,
    file.path(contam_dir, paste0("mean_contaminated_", p, ".rds"))
  )
}

clean_partition <- readRDS(file.path(clean_dir, "clean.rds"))
clean_tau_hat_list <- vector("list", length(clean_partition))

splittingrule <- "user"

colnames(lalonde_data) <- make.names(colnames(lalonde_data))

feature_cols_prop <- setdiff(colnames(lalonde_data), c("re78", "treat"))

propensity_model <- glm(
  as.formula(
    paste("treat ~", paste(feature_cols_prop, collapse = " + "))
  ),
  data = lalonde_data,
  family = binomial
)

for (i in 1:n_splits) {
  cat("Training clean model on split", i, "\n")
  
  datTrain <- clean_partition[[i]]$train
  datTest <- clean_partition[[i]]$test
  
  colnames(datTrain) <- make.names(colnames(datTrain))
  colnames(datTest) <- make.names(colnames(datTest))
  
  feature_cols <- setdiff(colnames(datTrain), c("re78", "treat"))
  f_2 <- paste(feature_cols, collapse = " + ")
  
  forest_clean <- causalTree::causalForest(
    as.formula(paste("re78 ~", f_2)),
    data = datTrain,
    treatment = datTrain$treat,
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
  clean_tau_hat_list[[i]] <- as.numeric(predict(forest_clean, Xtest))
}

saveRDS(
  clean_tau_hat_list,
  file.path(clean_results_dir, "tau_hat_clean_list.rds")
)

tau_hat_clean_list <- readRDS(
  file.path(clean_results_dir, "tau_hat_clean_list.rds")
)

for (p in p_values) {
  contaminated_splits <- readRDS(
    file.path(contam_dir, paste0("mean_contaminated_", p, ".rds"))
  )
  
  results_all <- vector("list", length(contaminated_splits))
  
  for (i in 1:n_splits) {
    cat("Processing contaminated split", i, "for p =", p, "\n")
    
    datTrain <- contaminated_splits[[i]]$train
    datTest <- contaminated_splits[[i]]$test
    
    colnames(datTrain) <- make.names(colnames(datTrain))
    colnames(datTest) <- make.names(colnames(datTest))
    
    feature_cols <- setdiff(colnames(datTrain), c("re78", "treat"))
    f_2 <- paste(feature_cols, collapse = " + ")
    
    Xtrain <- datTrain[, feature_cols]
    Xtest <- datTest[, feature_cols]
    
    forest_honest <- causalTree::causalForest(
      as.formula(paste("re78 ~", f_2)),
      data = datTrain,
      treatment = datTrain$treat,
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
      propensity_model, newdata = datTest, type = "response"
    )
    propens <- pmin(pmax(propens, 1e-6), 1 - 1e-6)
    
    ystar <- datTest$re78 *
      (datTest$treat - propens) / (propens * (1 - propens))
    
    tau_hat_clean <- tau_hat_clean_list[[i]]
    
    mape <- mean(
      abs((tau_hat - tau_hat_clean) / (abs(tau_hat_clean) + 1e-6))
    ) * 100
    
    rmse <- sqrt(mean((tau_hat - tau_hat_clean)^2, na.rm = TRUE))
    mae <- mean(abs(tau_hat - tau_hat_clean))
    cor_tau <- cor(tau_hat, tau_hat_clean)
    
    results_all[[i]] <- tibble(
      tau_hat = tau_hat,
      tau_hat_clean = tau_hat_clean,
      ystar = ystar,
      mape = mape * 100,
      rmse = rmse,
      mae = mae,
      correlation = cor_tau
    )
  }
  
  saveRDS(
    results_all,
    file.path(
      contam_results_dir,
      paste0("results_all_", p, ".rds")
    )
  )
}