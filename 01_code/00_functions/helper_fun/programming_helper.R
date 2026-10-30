# parallel function ----

number_of_workers <- function(n = NULL, ...) {
  if (!is.null(n) && !is.numeric(n)) {
    cli::cli_abort("Number of workers {.arg n} has to be {.cls numeric} or {.val NULL}.")
  }
  
  if (is.null(n)) {
    n <- floor(parallel::detectCores() * 0.5)
  } else if (n > parallel::detectCores()) {
    cli::cli_abort("Number of workers {.arg n} cannot be greater than the number of available CPU cores ({.val {parallel::detectCores()}}).")
  }
  
  cli::cli_alert_info("{.val {n}} workers set for parallel computing.")
  invisible(n)
}


load_rds_safely <- function(folder_path, file_name) {
  file_path <- here::here(folder_path, glue::glue("{file_name}.Rds"))
  
  if (!file.exists(file_path)) {
    cli::cli_abort(c(
      "File not found: {.file {file_name}.Rds}",
      "x" = "Expected location: {.path {file_path}}",
      "i" = "Please check that {.var folder_path} and {.var file_name} are correct."
    ))
  }
  
  readRDS(file_path)
}
