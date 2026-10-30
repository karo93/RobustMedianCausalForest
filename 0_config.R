###############################################################################
# @Org - University of Duisburg-Essen
# @Proj - RobustMedianCausalForest
# @Descrip - config file (declare glob vars, paths; load functions, ...)
###############################################################################


# Load here-package ----
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

# Load packages ----
source(here::here("00_packages/packages.R"))


# source all function  ----
invisible(
  sapply(
    list.files(
      here::here('01_code/00_functions'), recursive = TRUE, full.names = TRUE),
    source)
  )
