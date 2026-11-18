# Load config file ----
source(here::here("0_config.R"))

# parallel plan ----
# setting the number of workers (default 50%)
n_workers <- number_of_workers()
# paralletl plan
future::plan(multisession, workers = n_workers)
# Increase max allowed size for exported globals (2GB)
options(future.globals.maxSize = 2 * 1024^3)

# scenario ----
scenario <- 'no_effect'

#### Simulation settings ---
# Normal
dgp_setting_normal <- tibble::tibble(
  n_dataset = 10,  # number of datasets
  N = 2000, # number of obs. per sim. dataset
  D = 10, # number of covs per dataset
  noise = "normal", # error term noise
  noise_parameters = list(list('mean' = 0, 'sd' = 1)),
  tau_effect = "het", # treatment
  tau_variant = scenario,
  sigma_variant = 2, # sigmund variant
  x_dist = 'normal')

# Cauchy 
dgp_setting_cauchy <- dgp_setting_normal %>%
  dplyr::mutate(
    noise = "cauchy",
    noise_parameters = list(list('loc' = 0, 'scale' = 1))
    )

# t-distribution  
dgp_setting_heavy <- dgp_setting_normal %>%
  dplyr::mutate(
    noise = "heavy", 
    noise_parameters = list(list('df' = 3))
    )

# La Place
dgp_setting_laplace <- dgp_setting_normal %>%
  dplyr::mutate(
    noise = "laplace",
    noise_parameters = list(list('loc' = 0, 'scale' = 1))
    )

#### Running Simulation  ---

# Run Normal Simulation
run_sim(dgp_setting_normal)

# Run Cauchy Simulation
run_sim(dgp_setting_cauchy)

# Run t-distribution Simulation
run_sim(dgp_setting_heavy)

# Run La Place Simulation
run_sim(dgp_setting_laplace)

#### Analyzing Simulation per case ---

analysis_fun(dgp_setting_normal)

analysis_fun(dgp_setting_cauchy)

analysis_fun(dgp_setting_heavy)

analysis_fun(dgp_setting_laplace)


### Constructing Tables and Plots over all scenarios ----

combine_evaluation("03_analysis/01_sim")






# TABLE 1


# ci later BLB and IJ ---
# Table 4 & 5
# Plot 10 - 13


