
# Train the first causal forest with the clean data
cf_clean <- causalTree::causalForest(
  Y_train ~ X1 + X2, data = as.data.frame(data_train), treatment = W_train,
  split.Rule = "CT", cv.option = "CT", split.Honest = TRUE, cv.Honest = TRUE,
  split.Bucket = FALSE, ncov_sample = 2, ncolx = p, num.trees = 1000
)

# Train the second causal forest with a different splitting rule
cf_clean_alt <- causalTree::causalForest(
  Y_train ~ X1 + X2, data = as.data.frame(data_train), treatment = W_train,
  split.Rule = "user" , # Example alternative rule
  cv.option = "CT", split.Honest = TRUE, cv.Honest = TRUE,
  split.Bucket = FALSE, ncov_sample = 2, ncolx = p, num.trees = 1000
)

# Initialize variables to store errors
outlier <- relative_error_clean <- relative_error_clean_alt <- rep(NA, 10)
Y_outlier <- Y_train

# Loop through outliers
for (i in seq_along(outlier)) {
  Y_outlier[1] <- 100 * i + 2  # Introduce outlier
  data_train_outlier <- cbind(Y_outlier, X_train, W_train)
  
  # Train causal forest with outliers
  cf_outlier <- causalTree::causalForest(
    Y_outlier ~ X1 + X2, data = as.data.frame(data_train_outlier), treatment = W_train,
    split.Rule = "CT", cv.option = "CT", split.Honest = TRUE, cv.Honest = TRUE,
    split.Bucket = FALSE, ncov_sample = 2, ncolx = p, num.trees = 1000
  )
  
  cf_outlier_alt <- causalTree::causalForest(
    Y_outlier ~ X1 + X2, data = as.data.frame(data_train_outlier), treatment = W_train,
    split.Rule = "user", # Example alternative rule
    cv.option = "CT", split.Honest = TRUE, cv.Honest = TRUE,
    split.Bucket = FALSE, ncov_sample = 2, ncolx = p, num.trees = 1000
  )
  
  # Predict treatment effects
  tau_hat_clean <- predict(cf_clean, as.data.frame(X_test))
  tau_hat_clean_alt <- predict(cf_clean_alt, as.data.frame(X_test))
  tau_hat_outlier <- predict(cf_outlier, as.data.frame(X_test))
  tau_hat_outlier_alt <- predict(cf_outlier_alt, as.data.frame(X_test))
  
  # Compute relative prediction errors
  relative_error_clean[i] <- mean(abs(tau_hat_outlier[-1] - tau_test[-1])) / mean(abs(tau_hat_clean[-1] - tau_test[-1]))
  relative_error_clean_alt[i] <- mean(abs(tau_hat_outlier_alt[-1] - tau_test[-1])) / mean(abs(tau_hat_clean_alt[-1] - tau_test[-1]))
  outlier[i] <- Y_outlier[1]
}

# Combine the results in a data frame
df <- data.frame(
  OutlierValue = rep(outlier, 2),
  RelativeError = c(relative_error_clean, relative_error_clean_alt),
  Method = rep(c("CT Splitting Rule", "LMS Splitting Rule"), each = length(outlier))
)

# Plot the results
ggplot(df, aes(x = OutlierValue, y = RelativeError, color = Method)) +
  geom_line() +
  labs(
    y = "Mean prediction error",
    x = "Distance of the outlier"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("CT Splitting Rule" = "blue", "LMS Splitting Rule" = "red"))
