fit_single_pgs_model <- function(
    baseline,
    phenotype_row,
    pgs_row,
    covariates,
    current_use_var,
    population_override = NULL,
    analysis_family = "primary") {
  outcome_var <- phenotype_row$variable[[1]]
  pgs_var <- pgs_row$variable[[1]]
  model_family <- phenotype_row$family[[1]]
  population <- if (is.null(population_override)) {
    phenotype_row$primary_population[[1]]
  } else {
    population_override
  }

  required <- unique(c(outcome_var, pgs_var, covariates, current_use_var))
  analysis_data <- baseline[, ..required]
  if (population == "users_only") {
    analysis_data <- analysis_data[get(current_use_var) == 1]
  } else if (population != "all") {
    stop("Unknown primary_population: ", population, call. = FALSE)
  }
  analysis_data <- analysis_data[stats::complete.cases(analysis_data[, c(outcome_var, pgs_var, covariates), with = FALSE])]

  result_base <- data.table::data.table(
    analysis_family = analysis_family,
    phenotype = outcome_var,
    phenotype_label = phenotype_row$label[[1]],
    phenotype_domain = phenotype_row$domain[[1]],
    model_family = model_family,
    population = population,
    pgs = pgs_var,
    pgs_label = pgs_row$label[[1]],
    pgs_domain = pgs_row$domain[[1]],
    n = nrow(analysis_data)
  )

  if (nrow(analysis_data) < 50) {
    return(cbind(result_base, status = "insufficient_n", estimate = NA_real_, standard_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_))
  }

  analysis_data[, pgs_z := safe_z(get(pgs_var))]
  if (all(is.na(analysis_data$pgs_z))) {
    return(cbind(result_base, status = "constant_pgs", estimate = NA_real_, standard_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_))
  }
  analysis_data[, sex_model := factor(get(covariates[grepl("sex", covariates, ignore.case = TRUE)][1]))]
  model_covariates <- covariates
  sex_candidates <- model_covariates[grepl("sex", model_covariates, ignore.case = TRUE)]
  if (length(sex_candidates)) {
    model_covariates[model_covariates == sex_candidates[1]] <- "sex_model"
  }

  if (model_family == "continuous") {
    analysis_data[, outcome_model := safe_z(get(outcome_var))]
  } else if (model_family == "ordinal") {
    analysis_data[, outcome_model := ordered(get(outcome_var))]
  } else {
    analysis_data[, outcome_model := as.numeric(get(outcome_var))]
  }

  if (length(unique(stats::na.omit(analysis_data$outcome_model))) < 2) {
    return(cbind(result_base, status = "constant_outcome", estimate = NA_real_, standard_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_))
  }

  full_formula <- model_matrix_formula("outcome_model", c("pgs_z", model_covariates))
  reduced_formula <- model_matrix_formula("outcome_model", model_covariates)

  fit_warnings <- character()
  fitted <- tryCatch(withCallingHandlers({
    if (model_family == "continuous") {
      full <- stats::lm(full_formula, data = analysis_data)
      reduced <- stats::lm(reduced_formula, data = analysis_data)
      coefficient <- robust_coefficient(full, "pgs_z")
      fit_statistic <- summary(full)$adj.r.squared - summary(reduced)$adj.r.squared
      fit_statistic_name <- "incremental_adjusted_r2"
      effect_scale <- "SD outcome per 1-SD PGS"
      assumption_test <- NA_character_
      assumption_p_value <- NA_real_
      effect_ratio <- NA_real_
      ratio_low <- NA_real_
      ratio_high <- NA_real_
    } else if (model_family == "binary") {
      full <- stats::glm(full_formula, data = analysis_data, family = stats::binomial())
      reduced <- stats::glm(reduced_formula, data = analysis_data, family = stats::binomial())
      coefficient <- robust_coefficient(full, "pgs_z")
      fit_statistic <- tjur_r2(analysis_data$outcome_model, stats::fitted(full)) -
        tjur_r2(analysis_data$outcome_model, stats::fitted(reduced))
      fit_statistic_name <- "incremental_tjur_r2"
      effect_scale <- "Odds ratio per 1-SD PGS"
      assumption_test <- NA_character_
      assumption_p_value <- NA_real_
      effect_ratio <- exp(coefficient$estimate)
      ratio_low <- exp(coefficient$conf_low)
      ratio_high <- exp(coefficient$conf_high)
    } else if (model_family == "count") {
      full <- MASS::glm.nb(full_formula, data = analysis_data)
      reduced <- MASS::glm.nb(reduced_formula, data = analysis_data)
      nb_convergence_warning <- any(grepl(
        "iteration limit reached|NaNs produced|estimate truncated",
        fit_warnings
      ))
      if (nb_convergence_warning) {
        # The Poisson mean model with HC3 inference is a stable count-model
        # fallback when negative-binomial dispersion estimation is on its
        # numerical boundary. The warning is retained in the output audit.
        full <- stats::glm(
          full_formula, data = analysis_data, family = stats::poisson()
        )
        reduced <- stats::glm(
          reduced_formula, data = analysis_data, family = stats::poisson()
        )
      }
      coefficient <- robust_coefficient(full, "pgs_z")
      fit_statistic <- delta_loglik_metric(full, reduced)
      fit_statistic_name <- "incremental_loglik_fraction"
      effect_scale <- "Rate ratio per 1-SD PGS"
      assumption_test <- if (nb_convergence_warning) {
        "Poisson fallback with HC3 SE after negative-binomial convergence warning"
      } else {
        NA_character_
      }
      assumption_p_value <- NA_real_
      effect_ratio <- exp(coefficient$estimate)
      ratio_low <- exp(coefficient$conf_low)
      ratio_high <- exp(coefficient$conf_high)
    } else if (model_family == "ordinal") {
      full <- ordinal::clm(full_formula, data = analysis_data, link = "logit", Hess = TRUE)
      reduced <- ordinal::clm(reduced_formula, data = analysis_data, link = "logit", Hess = TRUE)
      coefficient_table <- coef(summary(full))
      if (!"pgs_z" %in% rownames(coefficient_table)) stop("PGS coefficient unavailable")
      estimate <- unname(coefficient_table["pgs_z", "Estimate"])
      standard_error <- unname(coefficient_table["pgs_z", "Std. Error"])
      p_column <- grep("Pr", colnames(coefficient_table), value = TRUE)[1]
      coefficient <- data.table::data.table(
        estimate = estimate,
        standard_error = standard_error,
        conf_low = estimate - stats::qnorm(0.975) * standard_error,
        conf_high = estimate + stats::qnorm(0.975) * standard_error,
        p_value = unname(coefficient_table["pgs_z", p_column])
      )
      fit_statistic <- delta_loglik_metric(full, reduced)
      fit_statistic_name <- "incremental_loglik_fraction"
      effect_scale <- "Proportional odds ratio per 1-SD PGS"
      nominal_test <- tryCatch(ordinal::nominal_test(full), error = function(error) NULL)
      assumption_test <- "PGS proportional-odds nominal test"
      assumption_p_value <- if (!is.null(nominal_test) && "pgs_z" %in% rownames(nominal_test)) {
        p_column <- grep("Pr", colnames(nominal_test), value = TRUE)[1]
        unname(nominal_test["pgs_z", p_column])
      } else {
        NA_real_
      }
      effect_ratio <- exp(coefficient$estimate)
      ratio_low <- exp(coefficient$conf_low)
      ratio_high <- exp(coefficient$conf_high)
    } else {
      stop("Unsupported model family: ", model_family)
    }
    if (is.null(coefficient)) stop("PGS coefficient unavailable")
    cbind(
      coefficient,
      effect_ratio = effect_ratio,
      ratio_low = ratio_low,
      ratio_high = ratio_high,
      effect_scale = effect_scale,
      fit_statistic_name = fit_statistic_name,
      fit_statistic = fit_statistic,
      assumption_test = assumption_test,
      assumption_p_value = assumption_p_value
    )
  }, warning = function(warning) {
    fit_warnings <<- c(fit_warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) {
    data.table::data.table(
      estimate = NA_real_, standard_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
      p_value = NA_real_, effect_ratio = NA_real_, ratio_low = NA_real_, ratio_high = NA_real_,
      effect_scale = NA_character_, fit_statistic_name = NA_character_, fit_statistic = NA_real_,
      assumption_test = NA_character_, assumption_p_value = NA_real_,
      error_message = conditionMessage(error)
    )
  })

  fitted$warning_message <- if (length(fit_warnings)) {
    paste(unique(fit_warnings), collapse = " | ")
  } else {
    NA_character_
  }

  status <- if (is.na(fitted$p_value[[1]])) "model_error" else "ok"
  cbind(result_base, status = status, fitted)
}

fit_pgs_manifest_screen <- function(
    baseline,
    pgs_manifest,
    phenotype_manifest,
    covariates,
    current_use_var,
    population_override = NULL,
    analysis_family) {
  primary_results <- vector("list", nrow(pgs_manifest) * nrow(phenotype_manifest))
  index <- 1L
  for (phenotype_index in seq_len(nrow(phenotype_manifest))) {
    phenotype_row <- phenotype_manifest[phenotype_index]
    for (pgs_index in seq_len(nrow(pgs_manifest))) {
      pgs_row <- pgs_manifest[pgs_index]
      primary_results[[index]] <- fit_single_pgs_model(
        baseline = baseline,
        phenotype_row = phenotype_row,
        pgs_row = pgs_row,
        covariates = covariates,
        current_use_var = current_use_var,
        population_override = population_override,
        analysis_family = analysis_family
      )
      index <- index + 1L
    }
  }
  data.table::rbindlist(primary_results, fill = TRUE)
}

build_global_hit_population_comparison <- function(
    baseline,
    primary,
    pgs_manifest,
    phenotype_manifest,
    covariates,
    current_use_var,
    excluded_phenotypes = c(
      "cannabis_ever", "cannabis_current",
      "logBreakpoint", "logIntensity", "logOmax", "logPmax", "logAlpha"
    )) {
  # Cannabis-use-status outcomes are constant among current users. MPT outcomes
  # are defined only among current users, so their nominal all-participant fits
  # use the same non-missing analysis sample and do not provide a meaningful
  # all-participants-versus-users comparison.
  selected <- primary[
    passes_global_fdr == TRUE & !phenotype %in% excluded_phenotypes
  ]
  if (!nrow(selected)) {
    empty <- data.table::copy(primary[0])
    empty[, `:=`(
      comparison_population = character(),
      estimate_role = character(),
      selection_p_value = numeric(),
      selection_fdr_q_global = numeric(),
      selected_by_global_fdr = logical(),
      inferential_status = character()
    )]
    return(empty)
  }

  data.table::rbindlist(lapply(seq_len(nrow(selected)), function(index) {
    primary_row <- data.table::copy(selected[index])
    phenotype_row <- phenotype_manifest[variable == primary_row$phenotype]
    pgs_row <- pgs_manifest[variable == primary_row$pgs]
    comparison_population <- if (primary_row$population == "all") "users_only" else "all"
    comparison_row <- fit_single_pgs_model(
      baseline, phenotype_row, pgs_row, covariates, current_use_var,
      population_override = comparison_population,
      analysis_family = "global_hit_population_effect_comparison"
    )
    primary_row[, `:=`(
      comparison_population = data.table::fifelse(
        population == "all", "All participants", "Current users only"
      ),
      estimate_role = "primary",
      inferential_status = "primary_global_fdr"
    )]
    comparison_row[, `:=`(
      comparison_population = data.table::fifelse(
        population == "all", "All participants", "Current users only"
      ),
      estimate_role = "effect_comparison",
      p_value = NA_real_,
      assumption_p_value = NA_real_,
      inferential_status = "effect_size_only_no_p_value"
    )]
    comparison <- data.table::rbindlist(list(primary_row, comparison_row), fill = TRUE)
    comparison[, `:=`(
      selection_p_value = primary_row$p_value,
      selection_fdr_q_global = primary_row$fdr_q_global,
      selected_by_global_fdr = TRUE
    )]
    comparison
  }), fill = TRUE)
}

run_baseline_pgs_screen <- function(
    baseline,
    pgs_manifest,
    phenotype_manifest,
    covariates,
    current_use_var,
    output_dir,
    global_fdr_alpha = 0.05) {
  primary <- fit_pgs_manifest_screen(
    baseline, pgs_manifest, phenotype_manifest, covariates, current_use_var,
    analysis_family = "primary_scale_appropriate"
  )
  primary[, fdr_q_global := stats::p.adjust(p_value, method = "BH", n = .N)]
  primary[, passes_global_fdr :=
    !is.na(fdr_q_global) & fdr_q_global < global_fdr_alpha]
  primary[, fdr_family_tests := .N]
  significant <- primary[passes_global_fdr == TRUE]

  # Structurally zero-coded outcomes are primary in all participants. Their
  # current-user versions are effect-size-only sensitivities and do not enter
  # the global FDR family.
  sensitivity_manifest <- phenotype_manifest[
    primary_population == "all" & zero_coded_sensitivity == TRUE
  ]
  user_sensitivity <- fit_pgs_manifest_screen(
    baseline, pgs_manifest, sensitivity_manifest, covariates, current_use_var,
    population_override = "users_only",
    analysis_family = "user_only_effect_sensitivity"
  )
  user_sensitivity[, `:=`(
    p_value = NA_real_,
    assumption_p_value = NA_real_,
    inferential_status = "effect_size_only_no_p_value"
  )]
  comparison <- build_global_hit_population_comparison(
    baseline, primary, pgs_manifest, phenotype_manifest, covariates,
    current_use_var
  )

  pgs_summary <- primary[, .(
    tests = .N,
    successful_models = sum(status == "ok"),
    fdr_significant_phenotypes = sum(passes_global_fdr, na.rm = TRUE),
    min_fit_statistic = suppressWarnings(min(fit_statistic[passes_global_fdr], na.rm = TRUE)),
    max_fit_statistic = suppressWarnings(max(fit_statistic[passes_global_fdr], na.rm = TRUE))
  ), by = .(pgs, pgs_label, pgs_domain)]
  pgs_summary[!is.finite(min_fit_statistic), min_fit_statistic := NA_real_]
  pgs_summary[!is.finite(max_fit_statistic), max_fit_statistic := NA_real_]

  phenotype_summary <- primary[, .(
    tests = .N,
    successful_models = sum(status == "ok"),
    fdr_significant_pgs = sum(passes_global_fdr, na.rm = TRUE)
  ), by = .(phenotype, phenotype_label, phenotype_domain, model_family, population)]

  write_csv(primary, file.path(output_dir, "baseline", "full_pgs_phenotype_results.csv"))
  write_csv(significant, file.path(output_dir, "baseline", "fdr_significant_results.csv"))
  write_csv(significant, file.path(output_dir, "baseline", "global_fdr_significant_results.csv"))
  write_csv(
    user_sensitivity,
    file.path(output_dir, "baseline", "user_only_effect_sensitivity.csv")
  )
  write_csv(comparison, file.path(output_dir, "baseline", "global_fdr_hits_all_vs_users_effects.csv"))
  write_csv(pgs_summary, file.path(output_dir, "baseline", "pgs_summary.csv"))
  write_csv(phenotype_summary, file.path(output_dir, "baseline", "phenotype_summary.csv"))

  list(
    primary = primary,
    significant = significant,
    user_sensitivity = user_sensitivity,
    population_comparison = comparison,
    sensitivity = user_sensitivity,
    pgs_summary = pgs_summary,
    phenotype_summary = phenotype_summary
  )
}
