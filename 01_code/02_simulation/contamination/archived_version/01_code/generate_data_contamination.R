data_contamination <- function(data, column, p, newValues) {
  numReplace <- ceiling(nrow(data) * p / 100)
  replaceIndices <- sample(1:nrow(data), size = numReplace)
  
  data[replaceIndices, column] <- sample(
    newValues, size = numReplace, replace = TRUE
  )
  
  return(data)
}

scenario <- "S4"

input_dir <- here::here("00_data", scenario, "00_Clean", "01_Train")
output_dir <- here::here("00_data", scenario, "01_Mean")

dataset_names <- list.files(input_dir)[1]
dataset_names_tools <- tools::file_path_sans_ext(dataset_names)

p_values <- c(5, 10, 15, 20, 25, 30, 35)

column_index <- "Yobs"
# sigma <- ...  # Set the value used for this scenario.

#table for mean/variance contamination
# S1: sigma=49.98,
# S2: sigma=1.15, 
# S3: sigma=17.23, 
# S4: sigma=8.45, 

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (p in p_values) {
  for (j in seq_along(dataset_names)) {
    
    dat1 <- readRDS(file.path(input_dir, dataset_names[j]))
    
    Ycontamin <- dat1$Yobs + rnorm(nrow(dat1), 5 * sigma, sigma)
    #Ycontamin <- dat1$Yobs + rnorm(nrow(dat1), 0, 5*sigma)  # Variance contamination

    dat1cont <- data_contamination(
      dat1,
      column = column_index,
      p = p,
      newValues = Ycontamin
    )
    
    saveRDS(
      dat1cont,
      file.path(
        output_dir,
        paste0(dataset_names_tools[j], "_", p, ".rds")
      )
    )
  }
}