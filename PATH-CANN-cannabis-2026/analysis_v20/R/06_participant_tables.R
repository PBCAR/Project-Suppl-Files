# PATH-CANN analysis v20 participant-table module.

normalize_table_sex <- function(x) {
  raw <- trimws(as.character(x))
  normalized <- rep(NA_character_, length(raw))
  lower <- tolower(raw)
  normalized[lower %in% c("male", "m", "man")] <- "Male"
  normalized[lower %in% c("female", "f", "woman")] <- "Female"

  numeric_value <- suppressWarnings(as.numeric(raw))
  numeric_observed <- sort(unique(numeric_value[!is.na(numeric_value)]))
  if (length(numeric_observed) && all(numeric_observed %in% c(1, 2))) {
    normalized[is.na(normalized) & numeric_value == 1] <- "Male"
    normalized[is.na(normalized) & numeric_value == 2] <- "Female"
  }

  unresolved <- !is.na(x) & is.na(normalized)
  if (any(unresolved)) {
    stop(
      "Sex-stratified tables require Male/Female labels or the audited 1=Male, 2=Female coding. Unresolved values: ",
      paste(unique(raw[unresolved]), collapse = ", "),
      call. = FALSE
    )
  }
  factor(normalized, levels = c("Female", "Male"))
}

table_population_data <- function(data, population, current_use_var) {
  if (population == "all") return(data)
  if (population == "current_users") {
    return(data[!is.na(get(current_use_var)) & get(current_use_var) == 1])
  }
  if (population == "ever_users") {
    return(data[!is.na(cannabis_ever) & cannabis_ever == 1])
  }
  stop("Unknown table population: ", population, call. = FALSE)
}

table_population_label <- function(population) {
  switch(
    population,
    all = "All participants",
    current_users = "Current cannabis users",
    ever_users = "Participants reporting lifetime cannabis use",
    stop("Unknown table population: ", population, call. = FALSE)
  )
}

format_table_number <- function(x, digits = 2L) {
  if (!is.finite(x)) return(NA_character_)
  formatC(x, digits = digits, format = "f", big.mark = ",")
}

summarize_table_vector <- function(x, summary_type, event_value = NA_character_) {
  numeric_x <- suppressWarnings(as.numeric(as.character(x)))
  observed <- numeric_x[!is.na(numeric_x)]
  n <- length(observed)
  missing <- length(x) - n
  if (!n) {
    return(list(
      statistic = NA_character_, n = 0L, missing = missing,
      observed_min = NA_real_, observed_max = NA_real_
    ))
  }

  statistic <- switch(
    summary_type,
    n = formatC(n, format = "d", big.mark = ","),
    n_percent = {
      event <- suppressWarnings(as.numeric(event_value))
      numerator <- sum(observed == event)
      paste0(
        formatC(numerator, format = "d", big.mark = ","),
        " (", format_table_number(100 * numerator / n, 1L), "%)"
      )
    },
    mean_sd = paste0(
      format_table_number(mean(observed), 2L), " (",
      format_table_number(stats::sd(observed), 2L), ")"
    ),
    median_iqr = {
      quantiles <- stats::quantile(observed, c(0.25, 0.5, 0.75), names = FALSE)
      paste0(
        format_table_number(quantiles[[2]], 2L), " [",
        format_table_number(quantiles[[1]], 2L), ", ",
        format_table_number(quantiles[[3]], 2L), "]"
      )
    },
    stop("Unsupported table summary type: ", summary_type, call. = FALSE)
  )
  list(
    statistic = statistic,
    n = n,
    missing = missing,
    observed_min = min(observed),
    observed_max = max(observed)
  )
}

summarize_table_category <- function(x, category_value) {
  observed <- as.character(x[!is.na(x)])
  n <- length(observed)
  missing <- length(x) - n
  numerator <- sum(observed == as.character(category_value))
  statistic <- if (n) {
    paste0(
      formatC(numerator, format = "d", big.mark = ","),
      " (", format_table_number(100 * numerator / n, 1L), "%)"
    )
  } else {
    NA_character_
  }
  list(statistic = statistic, n = n, missing = missing)
}

table_strata <- function(data) {
  list(
    Overall = data,
    Female = data[table_sex == "Female"],
    Male = data[table_sex == "Male"]
  )
}

category_definition <- function(variable, summary_type, data) {
  if (summary_type == "category" && variable == "ethnicity") {
    return(data.table::data.table(
      category_value = c("0", "1"),
      category_label = c("White", "Other racial group(s)")
    ))
  }
  values <- sort(unique(as.character(data[[variable]][!is.na(data[[variable]])])))
  data.table::data.table(category_value = values, category_label = values)
}

build_table_row <- function(
    data,
    section,
    variable,
    label,
    summary_type,
    population,
    event_value,
    scale,
    notes,
    current_use_var,
    family = NA_character_,
    domain = NA_character_,
    category_value = NA_character_,
    category_label = NA_character_) {
  population_data <- table_population_data(data, population, current_use_var)
  strata <- table_strata(population_data)
  row <- data.table::data.table(
    section = section,
    domain = domain,
    variable = variable,
    family = family,
    characteristic = if (!is.na(category_label)) {
      paste0(label, ": ", category_label)
    } else {
      label
    },
    scale_or_coding = scale,
    denominator_population = table_population_label(population),
    notes = notes
  )
  long_rows <- list()
  for (stratum_name in names(strata)) {
    stratum <- strata[[stratum_name]]
    summary <- if (variable == "__sample_n__") {
      list(
        statistic = formatC(nrow(stratum), format = "d", big.mark = ","),
        n = nrow(stratum), missing = 0L,
        observed_min = NA_real_, observed_max = NA_real_
      )
    } else if (summary_type %in% c("category", "category_auto")) {
      summarize_table_category(stratum[[variable]], category_value)
    } else {
      summarize_table_vector(stratum[[variable]], summary_type, event_value)
    }
    prefix <- tolower(stratum_name)
    row[, (paste0(prefix, "_statistic")) := summary$statistic]
    row[, (paste0(prefix, "_n")) := summary$n]
    row[, (paste0(prefix, "_missing")) := summary$missing]
    long_rows[[stratum_name]] <- data.table::data.table(
      section = section,
      domain = domain,
      variable = variable,
      family = family,
      characteristic = row$characteristic,
      scale_or_coding = scale,
      denominator_population = table_population_label(population),
      stratum = stratum_name,
      statistic = summary$statistic,
      n = summary$n,
      missing = summary$missing,
      observed_min = if (!is.null(summary$observed_min)) summary$observed_min else NA_real_,
      observed_max = if (!is.null(summary$observed_max)) summary$observed_max else NA_real_,
      notes = notes
    )
  }
  list(row = row, long = data.table::rbindlist(long_rows, fill = TRUE))
}

compare_table_sexes <- function(
    data,
    variable,
    summary_type,
    population,
    current_use_var) {
  population_data <- data.table::copy(
    table_population_data(data, population, current_use_var)
  )
  test_data <- population_data[
    !is.na(table_sex) & !is.na(get(variable)),
    .(
      table_sex,
      test_value = suppressWarnings(as.numeric(as.character(get(variable))))
    )
  ][!is.na(test_value)]
  female_n <- test_data[table_sex == "Female", .N]
  male_n <- test_data[table_sex == "Male", .N]
  result_base <- data.table::data.table(
    sex_difference_population = table_population_label(population),
    sex_difference_n = nrow(test_data),
    sex_difference_female_n = female_n,
    sex_difference_male_n = male_n
  )
  if (female_n == 0L || male_n == 0L) {
    return(cbind(
      result_base,
      sex_difference_test = NA_character_,
      sex_difference_statistic = NA_real_,
      sex_difference_df = NA_real_,
      sex_difference_p = NA_real_,
      sex_difference_status = "insufficient_sex_groups"
    ))
  }

  fitted <- tryCatch({
    if (summary_type == "mean_sd") {
      if (female_n < 2L || male_n < 2L) {
        stop("Welch test requires at least two observations in each sex group.")
      }
      test <- stats::t.test(
        test_value ~ table_sex, data = test_data, var.equal = FALSE
      )
      data.table::data.table(
        sex_difference_test = "Welch two-sample t-test",
        sex_difference_statistic = unname(test$statistic),
        sex_difference_df = unname(test$parameter),
        sex_difference_p = test$p.value
      )
    } else if (summary_type == "median_iqr") {
      test <- suppressWarnings(stats::wilcox.test(
        test_value ~ table_sex,
        data = test_data,
        exact = FALSE,
        correct = TRUE
      ))
      data.table::data.table(
        sex_difference_test = "Wilcoxon rank-sum test",
        sex_difference_statistic = unname(test$statistic),
        sex_difference_df = NA_real_,
        sex_difference_p = test$p.value
      )
    } else if (summary_type == "n_percent") {
      outcome_values <- sort(unique(test_data$test_value))
      if (length(outcome_values) != 2L) {
        stop("Binary sex comparison requires two observed outcome levels.")
      }
      contingency <- table(test_data$table_sex, test_data$test_value)
      chi <- suppressWarnings(stats::chisq.test(
        contingency, correct = FALSE
      ))
      if (any(chi$expected < 5)) {
        test <- stats::fisher.test(contingency)
        data.table::data.table(
          sex_difference_test = "Fisher exact test",
          sex_difference_statistic = NA_real_,
          sex_difference_df = NA_real_,
          sex_difference_p = test$p.value
        )
      } else {
        data.table::data.table(
          sex_difference_test = "Pearson chi-square test",
          sex_difference_statistic = unname(chi$statistic),
          sex_difference_df = unname(chi$parameter),
          sex_difference_p = chi$p.value
        )
      }
    } else {
      stop("No sex-comparison test is defined for summary type: ", summary_type)
    }
  }, error = function(error) {
    data.table::data.table(
      sex_difference_test = NA_character_,
      sex_difference_statistic = NA_real_,
      sex_difference_df = NA_real_,
      sex_difference_p = NA_real_,
      sex_difference_status = paste0("test_error: ", conditionMessage(error))
    )
  })
  if (!"sex_difference_status" %in% names(fitted)) {
    fitted[, sex_difference_status := data.table::fifelse(
      is.finite(sex_difference_p), "ok", "nonfinite_p_value"
    )]
  }
  cbind(result_base, fitted)
}

run_baseline_participant_tables <- function(
    baseline,
    phenotype_manifest,
    characteristics_manifest_path,
    current_use_var,
    sex_var,
    output_dir) {
  characteristics_manifest <- data.table::fread(characteristics_manifest_path)
  check_required_columns(
    characteristics_manifest,
    c(
      "order", "section", "variable", "label", "summary_type",
      "population", "event_value", "required", "scale", "notes"
    ),
    "Participant characteristics manifest"
  )
  data <- data.table::copy(baseline)
  data[, table_sex := normalize_table_sex(get(sex_var))]

  required_variables <- characteristics_manifest[
    required == TRUE & variable != "__sample_n__", variable
  ]
  check_required_columns(data, required_variables, "Baseline participant table data")
  ethnicity_value <- suppressWarnings(as.numeric(data$ethnicity))
  financial_value <- suppressWarnings(as.numeric(data$financial2))
  if (any(!is.na(ethnicity_value) & !ethnicity_value %in% c(0, 1))) {
    stop("ethnicity must retain the audited 0 = White, 1 = other race coding.", call. = FALSE)
  }
  if (any(!is.na(financial_value) & !financial_value %in% 1:9)) {
    stop("financial2 must be an ordered objective income bracket from 1 to 9.", call. = FALSE)
  }

  main_rows <- list()
  main_long <- list()
  row_index <- 1L
  for (manifest_index in seq_len(nrow(characteristics_manifest))) {
    spec <- characteristics_manifest[manifest_index]
    variable <- spec$variable[[1]]
    if (variable != "__sample_n__" && !variable %in% names(data)) next
    if (spec$summary_type[[1]] %in% c("category", "category_auto")) {
      categories <- category_definition(variable, spec$summary_type[[1]], data)
    } else {
      categories <- data.table::data.table(
        category_value = NA_character_, category_label = NA_character_
      )
    }
    for (category_index in seq_len(nrow(categories))) {
      built <- build_table_row(
        data,
        spec$section[[1]],
        variable,
        spec$label[[1]],
        spec$summary_type[[1]],
        spec$population[[1]],
        spec$event_value[[1]],
        spec$scale[[1]],
        spec$notes[[1]],
        current_use_var,
        category_value = categories$category_value[[category_index]],
        category_label = categories$category_label[[category_index]]
      )
      built$row[, display_order := spec$order[[1]] + category_index / 100]
      built$long[, display_order := spec$order[[1]] + category_index / 100]
      main_rows[[row_index]] <- built$row
      main_long[[row_index]] <- built$long
      row_index <- row_index + 1L
    }
  }
  main <- data.table::rbindlist(main_rows, fill = TRUE)
  main_long <- data.table::rbindlist(main_long, fill = TRUE)
  data.table::setorder(main, display_order)
  data.table::setorder(main_long, display_order, stratum)

  supplement_rows <- list()
  supplement_long <- list()
  for (index in seq_len(nrow(phenotype_manifest))) {
    spec <- phenotype_manifest[index]
    variable <- spec$variable[[1]]
    population <- if (variable == "cannabis_firstuse") {
      "ever_users"
    } else if (spec$primary_population[[1]] == "users_only") {
      "current_users"
    } else {
      spec$primary_population[[1]]
    }
    summary_type <- if (spec$family[[1]] == "binary") {
      "n_percent"
    } else if (
      spec$family[[1]] %in% c("count", "ordinal") ||
        variable == "cannabis_expenditure"
    ) {
      "median_iqr"
    } else {
      "mean_sd"
    }
    event_value <- if (summary_type == "n_percent") "1" else NA_character_
    scale <- spec$value_coding[[1]]
    legalization_note <- if (variable %in% c("nsduh3", "cannabis_risk")) {
      " Study item included to characterize legalization-context perceptions; scaling is reported explicitly."
    } else {
      ""
    }
    built <- build_table_row(
      data,
      spec$domain[[1]],
      variable,
      spec$label[[1]],
      summary_type,
      population,
      event_value,
      scale,
      paste0(spec$notes[[1]], legalization_note),
      current_use_var,
      family = spec$family[[1]],
      domain = spec$domain[[1]]
    )
    sex_test <- compare_table_sexes(
      data,
      variable,
      summary_type,
      population,
      current_use_var
    )
    for (test_column in names(sex_test)) {
      built$row[, (test_column) := sex_test[[test_column]][[1]]]
      built$long[, (test_column) := sex_test[[test_column]][[1]]]
    }
    built$row[, display_order := index]
    built$long[, display_order := index]
    supplement_rows[[index]] <- built$row
    supplement_long[[index]] <- built$long
  }
  supplement <- data.table::rbindlist(supplement_rows, fill = TRUE)
  supplement_long <- data.table::rbindlist(supplement_long, fill = TRUE)
  supplement[, sex_difference_fdr_q := stats::p.adjust(
    sex_difference_p, method = "BH", n = .N
  )]
  supplement_long <- merge(
    supplement_long,
    supplement[, .(variable, sex_difference_fdr_q)],
    by = "variable", all.x = TRUE, sort = FALSE
  )

  table_notes <- data.table::data.table(
    note_order = 1:13,
    note = c(
      "Values are mean (SD), median [Q1, Q3], or n (%), as indicated by the characteristic and scale.",
      "Every row reports its variable-specific non-missing denominator.",
      "Sex-stratified columns use the audited source coding 1 = Male and 2 = Female, or explicit Male/Female labels.",
      "Sex comparisons use exactly the population in which each phenotype is defined: all participants, current cannabis users, or participants reporting lifetime cannabis use.",
      "Welch two-sample t-tests are used for variables summarized by mean (SD), Wilcoxon rank-sum tests for variables summarized by median [Q1, Q3], and Pearson chi-square tests for binary variables; Fisher exact tests replace chi-square tests when any expected cell is below five.",
      "Raw sex-comparison p-values and Benjamini-Hochberg q-values across the 30 phenotype comparisons are reported; these comparisons are descriptive and exploratory.",
      "The variable historically named ethnicity is an audited White-versus-other race indicator and is reported as race, not ethnicity.",
      "The ten ancestry principal components are adjustment covariates but are not interpretable participant characteristics and are therefore not tabulated individually.",
      "User-only phenotypes are summarized among current users; eligible-user missing responses remain missing and are never converted to zero.",
      "CUDIT-C, CUDIT-P, MMM motives, and expenditure use the analysis-ready structural-zero rule for confirmed non-users in all-participant summaries.",
      "Perceived weekly harm is coded 1 = no risk to 4 = great risk; perceived cannabis risk is a 0–12 count of endorsed risk options. Both characterize legalization-context perceptions.",
      "For medical and recreational cannabis acceptability, source code 5 = No opinion is treated as missing; the analyzed ordinal scale is 1–4.",
      "MPT parameters and cannabis-effect measures are defined among current users; non-users are excluded from these analyses. MPT parameters are recomputed as log1p of their raw non-negative values."
    )
  )

  write_csv(
    main,
    file.path(output_dir, "tables", "supplementary_table_s1_analytic_sample_characteristics.csv")
  )
  write_csv(
    main_long,
    file.path(output_dir, "tables", "supplementary_table_s1_analytic_sample_characteristics_long.csv")
  )
  write_csv(
    supplement,
    file.path(output_dir, "tables", "supplementary_table_s2_all_cannabis_phenotypes.csv")
  )
  write_csv(
    supplement_long,
    file.path(output_dir, "tables", "supplementary_table_s2_all_cannabis_phenotypes_long.csv")
  )
  write_csv(table_notes, file.path(output_dir, "tables", "baseline_table_notes.csv"))
  write_display_text(
    file.path(output_dir, "tables"),
    "supplementary_table_s1_analytic_sample_characteristics",
    "Table S1. Baseline characteristics of the analytic sample.",
    paste0(
      "Values are mean (SD), median [Q1, Q3], or n (%), with variable-specific ",
      "non-missing denominators. Combined and sex-stratified summaries are descriptive; ",
      "no sex-comparison p-values or standardized mean differences are reported."
    )
  )
  write_display_text(
    file.path(output_dir, "tables"),
    "supplementary_table_s2_all_cannabis_phenotypes",
    "Table S2. Baseline distributions of all cannabis-related phenotypes.",
    paste0(
      "Values are mean (SD), median [Q1, Q3], or n (%), with variable-specific ",
      "denominators in each phenotype's designated analysis population. Combined, Female ",
      "and Male summaries are accompanied by an appropriate two-sample sex comparison. ",
      "Welch t-tests are used for mean/SD variables, Wilcoxon rank-sum tests for median/IQR ",
      "variables, and Pearson chi-square or Fisher exact tests for binary variables. Raw p-values ",
      "and BH-adjusted q-values across all 30 comparisons are reported. Purpose-built legalization-",
      "context perception items retain their explicit scoring definitions."
    )
  )

  list(
    main = main,
    main_long = main_long,
    supplement = supplement,
    supplement_long = supplement_long,
    notes = table_notes
  )
}
