# Load CRAN-packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load("partykit",
               "pak",
               "data.tree",
               "Matching",
               "dplyr", 
               "devtools", 
               "randomForest",
               "rpart", 
               "grf", 
               "ggplot2", 
               "SimDesign", 
               "here", 
               "glue", 
               "future",
               "tidyr",
               "magrittr",
               "rpart.plot",
               "tictoc",
               "ranger",
               "ExtDist",
               "crayon",
               "cli"
)


# robustCausalTree-package
# changed for our purposes
if (!require("robustCausalTree", quietly = TRUE)){
  # install local version
  install.packages(here::here("00_packages/robustCausalTree/"),
                   repos = NULL,
                   type = "source",
                   clean = TRUE)
}

library(robustCausalTree)