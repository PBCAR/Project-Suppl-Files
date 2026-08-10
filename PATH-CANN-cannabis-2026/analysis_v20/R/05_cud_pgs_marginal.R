# PATH-CANN analysis v20 marginal repeated-CUD module.
#
# Fits only the selected-PGS, wave-saturated marginal logistic models needed
# for the revised Figure 3 and its omnibus PGS-by-wave heterogeneity tests.
# No random effects are estimated.

marginal_interaction_term_name <- function(coefficient_names, wave_level) {
  candidates <- c(
    paste0("pgs_z:cud_wave", wave_level),
    paste0("cud_wave", wave_level, ":pgs_z")
  )
  found <- intersect(candidates, coefficient_names)
  if (length(found)) found[[1]] else NA_character_
}

prepare_marginal_cud_data <- function(
    baseline,
    longitudinal,
    pgs_var,
    waves,
    id_var,
    sex_var,
    wave_var,
    wave_age_var,
    cud_var,
    pc_vars) {
  optional_covariates <- intersect("genotyping_batch", names(baseline))
  baseline_columns <- unique(c(
    id_var, sex_var, pc_vars, optional_covariates, pgs_var
  ))
  check_required_columns(
    baseline, baseline_columns, paste0("Baseline data for ", pgs_var)
  )
  followup_columns <- c(id_var, wave_var, wave_age_var, cud_var)
  check_required_columns(
    longitudinal, followup_columns, paste0("Longitudinal data for ", pgs_var)
  )

  baseline_data <- data.table::copy(baseline[, ..baseline_columns])
  baseline_data[, `:=`(
    pgs_z = safe_z(get(pgs_var)),
    sex_model = factor(get(sex_var))
  )]
  if (length(optional_covariates)) {
    baseline_data[, genotyping_batch_model := factor(genotyping_batch)]
  }

  followup <- data.table::copy(longitudinal[
    get(wave_var) %in% waves,
    ..followup_columns
  ])
  data.table::setnames(
    followup,
    c(wave_var, wave_age_var, cud_var),
    c("cud_wave", "age_at_wave", "cud_status")
  )
  data <- merge(baseline_data, followup, by = id_var, all = FALSE)
  data[, `:=`(
    cud_wave = factor(cud_wave, levels = waves),
    cud_status = as.numeric(cud_status)
  )]
  required <- c(
    id_var, "cud_status", "cud_wave", "age_at_wave", "sex_model",
    "pgs_z", pc_vars
  )
  if (length(optional_covariates)) {
    required <- c(required, "genotyping_batch_model")
  }
  data <- data[stats::complete.cases(data[, ..required])]
  data[, age_z := safe_z(age_at_wave), by = cud_wave]
  data <- data[!is.na(age_z)]
  assert_binary(data$cud_status, "Repeated DSM-5 CUD status")
  data
}

extract_marginal_wave_effects <- function(
    model, covariance, pgs_row, waves, data, id_var) {
  beta <- stats::coef(model)
  coefficient_names <- names(beta)
  data.table::rbindlist(lapply(waves, function(wave) {
    contrast <- rep(0, length(beta))
    names(contrast) <- coefficient_names
    contrast[["pgs_z"]] <- 1
    if (wave != waves[[1]]) {
      interaction <- marginal_interaction_term_name(coefficient_names, wave)
      if (is.na(interaction)) {
        stop("Marginal PGS-by-wave coefficient unavailable for ", wave, call. = FALSE)
      }
      contrast[[interaction]] <- 1
    }
    estimate <- sum(contrast * beta)
    standard_error <- sqrt(drop(t(contrast) %*% covariance %*% contrast))
    statistic <- estimate / standard_error
    wave_data <- data[cud_wave == wave]
    data.table::data.table(
      pgs = pgs_row$pgs[[1]],
      pgs_label = pgs_row$pgs_label[[1]],
      pgs_domain = pgs_row$pgs_domain[[1]],
      wave = wave,
      n = nrow(wave_data),
      cases = sum(wave_data$cud_status == 1),
      noncases = sum(wave_data$cud_status == 0),
      estimate = estimate,
      standard_error = standard_error,
      conf_low = estimate - stats::qnorm(0.975) * standard_error,
      conf_high = estimate + stats::qnorm(0.975) * standard_error,
      p_value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
      effect_ratio = exp(estimate),
      ratio_low = exp(estimate - stats::qnorm(0.975) * standard_error),
      ratio_high = exp(estimate + stats::qnorm(0.975) * standard_error),
      effect_scale = "Marginal odds ratio per 1-SD PGS",
      n_participants_model = data.table::uniqueN(data[[id_var]]),
      participant_waves_model = nrow(data)
    )
  }))
}

extract_marginal_interaction_test <- function(
    model, covariance, pgs_row, waves, data, id_var) {
  coefficient_names <- names(stats::coef(model))
  interaction_terms <- vapply(
    waves[-1],
    function(wave) marginal_interaction_term_name(coefficient_names, wave),
    character(1)
  )
  if (anyNA(interaction_terms)) {
    stop(
      "Marginal PGS-by-wave interaction terms unavailable for: ",
      paste(waves[-1][is.na(interaction_terms)], collapse = ", "),
      call. = FALSE
    )
  }
  interaction_beta <- stats::coef(model)[interaction_terms]
  interaction_covariance <- covariance[
    interaction_terms, interaction_terms, drop = FALSE
  ]
  inverse_covariance <- tryCatch(
    solve(interaction_covariance),
    error = function(error) qr.solve(interaction_covariance)
  )
  statistic <- drop(
    t(interaction_beta) %*% inverse_covariance %*% interaction_beta
  )
  degrees_freedom <- length(interaction_terms)
  data.table::data.table(
    pgs = pgs_row$pgs[[1]],
    pgs_label = pgs_row$pgs_label[[1]],
    pgs_domain = pgs_row$pgs_domain[[1]],
    n_participants = data.table::uniqueN(data[[id_var]]),
    participant_waves = nrow(data),
    cud_positive_observations = sum(data$cud_status == 1),
    chi_square = statistic,
    degrees_freedom = degrees_freedom,
    p_value = stats::pchisq(statistic, degrees_freedom, lower.tail = FALSE)
  )
}

fit_selected_cud_pgs_marginal <- function(
    baseline,
    longitudinal,
    pgs_row,
    waves,
    id_var,
    sex_var,
    wave_var,
    wave_age_var,
    cud_var,
    pc_vars) {
  pgs_var <- pgs_row$pgs[[1]]
  data <- prepare_marginal_cud_data(
    baseline, longitudinal, pgs_var, waves, id_var, sex_var,
    wave_var, wave_age_var, cud_var, pc_vars
  )
  optional_covariates <- if (
    "genotyping_batch_model" %in% names(data)
  ) "genotyping_batch_model" else character()
  stratified_covariates <- c(
    "pgs_z", "age_z", "sex_model", pc_vars, optional_covariates
  )
  formula <- stats::as.formula(paste0(
    "cud_status ~ cud_wave * (",
    paste(stratified_covariates, collapse = " + "),
    ")"
  ))

  captured <- tryCatch({
    model <- stats::glm(formula, data = data, family = stats::binomial())
    covariance <- sandwich::vcovCL(
      model, cluster = data[[id_var]], type = "HC1"
    )
    list(
      wave_effects = extract_marginal_wave_effects(
        model, covariance, pgs_row, waves, data, id_var
      ),
      interaction = extract_marginal_interaction_test(
        model, covariance, pgs_row, waves, data, id_var
      ),
      formula = paste(deparse(formula), collapse = " ")
    )
  }, error = function(error) error)

  if (inherits(captured, "error")) {
    error_row <- data.table::data.table(
      pgs = pgs_row$pgs[[1]],
      pgs_label = pgs_row$pgs_label[[1]],
      pgs_domain = pgs_row$pgs_domain[[1]],
      status = "model_error",
      error_message = conditionMessage(captured)
    )
    return(list(wave_effects = error_row, interaction = error_row))
  }
  captured$wave_effects[, `:=`(
    status = "ok",
    model_formula = captured$formula,
    error_message = NA_character_
  )]
  captured$interaction[, `:=`(
    status = "ok",
    model_formula = captured$formula,
    error_message = NA_character_
  )]
  list(
    wave_effects = captured$wave_effects,
    interaction = captured$interaction
  )
}

run_selected_cud_pgs_marginal <- function(
    baseline,
    longitudinal,
    selection_audit,
    waves,
    id_var,
    sex_var,
    wave_var,
    wave_age_var,
    cud_var,
    pc_vars,
    output_dir,
    fdr_alpha = 0.05) {
  check_required_columns(
    selection_audit,
    c(
      "pgs", "pgs_label", "pgs_domain", "selected_for_cud_followup",
      "selection_reason"
    ),
    "Selected-PGS audit"
  )
  selected <- data.table::copy(
    selection_audit[selected_for_cud_followup %in% TRUE]
  )
  if (!nrow(selected)) {
    stop("No PGSs were selected for the marginal CUD figure.", call. = FALSE)
  }

  fitted <- lapply(seq_len(nrow(selected)), function(index) {
    fit_selected_cud_pgs_marginal(
      baseline, longitudinal, selected[index], waves, id_var, sex_var,
      wave_var, wave_age_var, cud_var, pc_vars
    )
  })
  wave_effects <- data.table::rbindlist(
    lapply(fitted, `[[`, "wave_effects"), fill = TRUE
  )
  interactions <- data.table::rbindlist(
    lapply(fitted, `[[`, "interaction"), fill = TRUE
  )
  wave_effects[status == "ok", `:=`(
    fdr_q_wave = stats::p.adjust(p_value, method = "BH", n = .N),
    fdr_family_tests = .N,
    inferential_status = "exploratory_selected_pgs_wave_specific_marginal"
  )]
  wave_effects[, passes_wave_fdr :=
    !is.na(fdr_q_wave) & fdr_q_wave < fdr_alpha]
  interactions[status == "ok", `:=`(
    fdr_q_interaction = stats::p.adjust(p_value, method = "BH", n = .N),
    fdr_family_tests = .N,
    inferential_status = "exploratory_selected_pgs_wave_interaction_marginal"
  )]
  interactions[, passes_interaction_fdr :=
    !is.na(fdr_q_interaction) & fdr_q_interaction < fdr_alpha]

  marginal_dir <- file.path(output_dir, "longitudinal")
  dir.create(marginal_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(
    wave_effects,
    file.path(marginal_dir, "cud_selected_pgs_wave_effects_marginal.csv")
  )
  write_csv(
    interactions,
    file.path(marginal_dir, "cud_selected_pgs_wave_interactions_marginal.csv")
  )
  list(wave_effects = wave_effects, interactions = interactions)
}
