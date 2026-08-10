required_packages <- c(
  "data.table", "broom", "sandwich", "lmtest", "MASS", "ordinal",
  "ggplot2", "patchwork", "scales"
)

check_packages <- function(packages = required_packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "),
      ". Install them before running this workflow.",
      call. = FALSE
    )
  }
}

read_analysis_file <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
  extension <- tolower(tools::file_ext(path))
  if (extension == "rds") {
    value <- readRDS(path)
  } else if (extension == "csv") {
    value <- data.table::fread(path)
  } else {
    stop("Unsupported input type for ", path, ". Use .rds or .csv.", call. = FALSE)
  }
  data.table::as.data.table(value)
}

safe_z <- function(x) {
  x <- as.numeric(x)
  value_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(value_sd) || value_sd == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / value_sd
}

clamp_probability <- function(x, epsilon = 1e-6) {
  pmin(pmax(as.numeric(x), epsilon), 1 - epsilon)
}

assert_unique_ids <- function(data, id_column, context) {
  duplicates <- data[!is.na(get(id_column)), .N, by = id_column][N > 1]
  if (nrow(duplicates)) {
    stop(context, " must have one row per participant; duplicated IDs found.", call. = FALSE)
  }
}

assert_binary <- function(x, name) {
  values <- sort(unique(stats::na.omit(x)))
  if (!all(values %in% c(0, 1))) {
    stop(name, " must be coded 0/1/NA; observed values: ", paste(values, collapse = ", "), call. = FALSE)
  }
}

check_required_columns <- function(data, required, context) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      context, " is missing required columns:\n- ",
      paste(missing, collapse = "\n- "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

initialize_output <- function(output_dir, run_id, config_snapshot) {
  if (dir.exists(output_dir)) {
    stop(
      "Output run already exists and will not be overwritten: ", output_dir,
      "\nChoose a new PATH_CANN_RUN_ID.",
      call. = FALSE
    )
  }
  subdirs <- c(
    "baseline", "longitudinal", "liability_bridge", "descriptives",
    "missingness", "figures", "tables", "qc"
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (subdir in subdirs) {
    dir.create(file.path(output_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  }
  log_lines <- c(
    paste0("Run ID: ", run_id),
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "Workflow version: analysis_v20",
    "Outputs are additive and this directory will not be overwritten.",
    "",
    "Configuration:",
    config_snapshot,
    "",
    paste0("R version: ", R.version.string)
  )
  writeLines(log_lines, file.path(output_dir, "version.log"), useBytes = TRUE)
}

write_csv <- function(data, path) {
  data.table::fwrite(data.table::as.data.table(data), path, na = "")
  invisible(path)
}

model_matrix_formula <- function(outcome, predictors) {
  stats::as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
}

robust_coefficient <- function(model, term, type = "HC3", cluster = NULL) {
  vc <- if (is.null(cluster)) {
    sandwich::vcovHC(model, type = type)
  } else {
    sandwich::vcovCL(model, cluster = cluster, type = "HC1")
  }
  table <- lmtest::coeftest(model, vcov. = vc)
  if (!term %in% rownames(table)) {
    return(NULL)
  }
  estimate <- unname(table[term, 1])
  standard_error <- unname(table[term, 2])
  p_value <- unname(table[term, 4])
  data.table::data.table(
    estimate = estimate,
    standard_error = standard_error,
    conf_low = estimate - stats::qnorm(0.975) * standard_error,
    conf_high = estimate + stats::qnorm(0.975) * standard_error,
    p_value = p_value
  )
}

tjur_r2 <- function(y, fitted) {
  if (length(unique(stats::na.omit(y))) < 2) return(NA_real_)
  mean(fitted[y == 1], na.rm = TRUE) - mean(fitted[y == 0], na.rm = TRUE)
}

delta_loglik_metric <- function(full_model, reduced_model) {
  full_ll <- as.numeric(stats::logLik(full_model))
  reduced_ll <- as.numeric(stats::logLik(reduced_model))
  if (!is.finite(full_ll) || !is.finite(reduced_ll) || reduced_ll == 0) return(NA_real_)
  (full_ll - reduced_ll) / abs(reduced_ll)
}

standardized_mean_difference <- function(x, included) {
  included <- as.logical(included)
  if (is.numeric(x) || is.integer(x)) {
    group_1 <- as.numeric(x[included])
    group_0 <- as.numeric(x[!included])
    pooled <- sqrt((stats::var(group_1, na.rm = TRUE) + stats::var(group_0, na.rm = TRUE)) / 2)
    if (!is.finite(pooled) || pooled == 0) return(NA_real_)
    return((mean(group_1, na.rm = TRUE) - mean(group_0, na.rm = TRUE)) / pooled)
  }
  x_factor <- factor(x)
  if (nlevels(x_factor) != 2) return(NA_real_)
  numeric_x <- as.integer(x_factor == levels(x_factor)[2])
  p1 <- mean(numeric_x[included], na.rm = TRUE)
  p0 <- mean(numeric_x[!included], na.rm = TRUE)
  denominator <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  (p1 - p0) / denominator
}

wave_number <- function(x) {
  suppressWarnings(as.integer(sub("^[^0-9]+", "", as.character(x))))
}
