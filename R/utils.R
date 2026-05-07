# check required column and optionally filter by minimum value
check_cols <- function(data, colname, min_value = NULL) {
  if (is.null(data)) return(NULL)
  if (!(colname %in% names(data))) return(NULL)
  if (!is.null(min_value)) {
    data <- data[data[[colname]] >= min_value, ]
    if (!nrow(data)) return(NULL)
  }
  data
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}