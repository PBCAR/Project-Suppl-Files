apply_v20_baseline_coding <- function(baseline, output_dir) {
  data <- data.table::copy(baseline)

  acceptability_variables <- c(
    "cann_med_acceptable", "cann_rec_acceptable"
  )
  check_required_columns(
    data, acceptability_variables,
    "Baseline data for acceptability recoding"
  )
  acceptability_audit <- data.table::rbindlist(lapply(
    acceptability_variables,
    function(variable) {
      source_value <- suppressWarnings(as.numeric(data[[variable]]))
      conversion_failure <- !is.na(data[[variable]]) & is.na(source_value)
      invalid <- !is.na(source_value) & !source_value %in% 1:5
      if (any(conversion_failure) || any(invalid)) {
        stop(
          variable,
          " must use source codes 1-5/NA before no-opinion recoding.",
          call. = FALSE
        )
      }
      no_opinion <- !is.na(source_value) & source_value == 5
      source_value[no_opinion] <- NA_real_
      data[[variable]] <<- source_value
      data.table::data.table(
        variable = variable,
        source_nonmissing = sum(!is.na(baseline[[variable]])),
        no_opinion_recode_to_missing = sum(no_opinion),
        analysis_nonmissing = sum(!is.na(source_value)),
        analysis_min = if (all(is.na(source_value))) NA_real_ else min(source_value, na.rm = TRUE),
        analysis_max = if (all(is.na(source_value))) NA_real_ else max(source_value, na.rm = TRUE)
      )
    }
  ))

  mpt_map <- data.table::data.table(
    raw_variable = c("Breakpoint", "Intensity", "Omax", "Pmax", "Alpha"),
    analysis_variable = c(
      "logBreakpoint", "logIntensity", "logOmax", "logPmax", "logAlpha"
    )
  )
  check_required_columns(
    data, mpt_map$raw_variable,
    "Baseline data for MPT log1p transformations"
  )
  mpt_audit <- data.table::rbindlist(lapply(
    seq_len(nrow(mpt_map)),
    function(index) {
      raw_variable <- mpt_map$raw_variable[[index]]
      analysis_variable <- mpt_map$analysis_variable[[index]]
      raw_source <- data[[raw_variable]]
      raw_value <- suppressWarnings(as.numeric(raw_source))
      conversion_failure <- !is.na(raw_source) & is.na(raw_value)
      if (any(conversion_failure)) {
        stop(raw_variable, " contains non-numeric MPT values.", call. = FALSE)
      }
      if (any(raw_value < 0, na.rm = TRUE)) {
        stop(
          raw_variable,
          " contains negative values and cannot use the prespecified log1p rule.",
          call. = FALSE
        )
      }
      transformed <- log1p(raw_value)
      data[[analysis_variable]] <<- transformed
      data.table::data.table(
        raw_variable = raw_variable,
        analysis_variable = analysis_variable,
        rule = paste0(analysis_variable, " = log1p(", raw_variable, ")"),
        raw_nonmissing = sum(!is.na(raw_value)),
        raw_zero = sum(raw_value == 0, na.rm = TRUE),
        transformed_nonmissing = sum(!is.na(transformed)),
        transformed_min = if (all(is.na(transformed))) NA_real_ else min(transformed, na.rm = TRUE),
        transformed_max = if (all(is.na(transformed))) NA_real_ else max(transformed, na.rm = TRUE)
      )
    }
  ))

  write_csv(
    acceptability_audit,
    file.path(output_dir, "qc", "acceptability_no_opinion_recode_audit.csv")
  )
  write_csv(
    mpt_audit,
    file.path(output_dir, "qc", "mpt_log1p_transformation_audit.csv")
  )
  list(
    data = data,
    acceptability_audit = acceptability_audit,
    mpt_audit = mpt_audit
  )
}

load_and_validate_inputs <- function(
    baseline_file,
    baseline_longitudinal_file,
    longitudinal_file,
    pgs_manifest_path,
    phenotype_manifest_path,
    id_var,
    baseline_age_var,
    wave_age_var,
    sex_var,
    wave_var,
    cud_var,
    cud_symptom_var,
    cud_symptom_observed_var,
    current_use_var,
    pc_vars,
    output_dir) {
  baseline <- read_analysis_file(baseline_file)
  baseline_longitudinal <- read_analysis_file(baseline_longitudinal_file)
  longitudinal <- read_analysis_file(longitudinal_file)
  pgs_manifest <- data.table::fread(pgs_manifest_path)
  phenotype_manifest <- data.table::fread(phenotype_manifest_path)

  check_required_columns(pgs_manifest, c("variable", "label", "domain"), "PGS manifest")
  check_required_columns(
    phenotype_manifest,
    c("variable", "label", "domain", "family", "primary_population", "zero_coded_sensitivity"),
    "Phenotype manifest"
  )

  coding <- apply_v20_baseline_coding(baseline, output_dir)
  baseline <- coding$data

  baseline_required <- unique(c(
    id_var, baseline_age_var, sex_var, pc_vars, current_use_var,
    "cannabis_frequency",
    pgs_manifest$variable, phenotype_manifest$variable
  ))
  longitudinal_required <- unique(c(
    id_var, wave_var, cud_var, cud_symptom_var, cud_symptom_observed_var,
    wave_age_var, sex_var,
    current_use_var, "cannabis_frequency", "cudit_freq", "cudit_sev",
    "mmm_coping", "nsduh3"
  ))
  baseline_longitudinal_required <- unique(c(
    id_var, baseline_age_var, sex_var, current_use_var,
    "cannabis_frequency", "cudit_freq", "cudit_sev", "mmm_coping", "nsduh3"
  ))

  check_required_columns(baseline, baseline_required, "Baseline analysis dataset")
  check_required_columns(
    baseline_longitudinal,
    baseline_longitudinal_required,
    "Full phenotype baseline predictor dataset"
  )
  check_required_columns(longitudinal, longitudinal_required, "Longitudinal analysis dataset")
  assert_unique_ids(baseline, id_var, "Baseline analysis dataset")
  assert_unique_ids(
    baseline_longitudinal,
    id_var,
    "Full phenotype baseline predictor dataset"
  )

  duplicate_long <- longitudinal[, .N, by = c(id_var, wave_var)][N > 1]
  if (nrow(duplicate_long)) {
    stop("Longitudinal data must have at most one row per participant-wave.", call. = FALSE)
  }
  assert_binary(baseline[[current_use_var]], current_use_var)
  assert_binary(
    baseline_longitudinal[[current_use_var]],
    paste0("full phenotype baseline ", current_use_var)
  )
  assert_binary(longitudinal[[current_use_var]], paste0("longitudinal ", current_use_var))
  assert_binary(longitudinal[[cud_var]], cud_var)
  symptom_count <- suppressWarnings(as.numeric(longitudinal[[cud_symptom_var]]))
  symptom_observed <- suppressWarnings(as.numeric(
    longitudinal[[cud_symptom_observed_var]]
  ))
  invalid_symptom_count <- !is.na(symptom_count) &
    (symptom_count < 0 | symptom_count > 11 | symptom_count %% 1 != 0)
  invalid_symptom_observed <- !is.na(symptom_observed) &
    (symptom_observed < 0 | symptom_observed > 11 | symptom_observed %% 1 != 0)
  invalid_partial_sum <- !is.na(symptom_count) & !is.na(symptom_observed) &
    symptom_count > symptom_observed
  if (any(invalid_symptom_count)) {
    stop(cud_symptom_var, " must be an integer from 0 to 11 or NA.", call. = FALSE)
  }
  if (any(invalid_symptom_observed)) {
    stop(
      cud_symptom_observed_var,
      " must be an integer from 0 to 11 or NA.",
      call. = FALSE
    )
  }
  if (any(invalid_partial_sum)) {
    stop(
      cud_symptom_var, " cannot exceed ", cud_symptom_observed_var, ".",
      call. = FALSE
    )
  }

  frequency <- suppressWarnings(as.numeric(baseline[["cannabis_frequency"]]))
  current <- suppressWarnings(as.numeric(baseline[[current_use_var]]))
  invalid_frequency <- !is.na(frequency) & !(frequency %in% 0:7)
  if (any(invalid_frequency)) {
    stop("cannabis_frequency must be coded 0-7/NA.", call. = FALSE)
  }
  eligibility_observed <- !is.na(frequency) & !is.na(current)
  eligibility_conflict <- eligibility_observed & (
    (current == 0 & frequency != 0) |
      (current == 1 & frequency == 0)
  )
  if (any(eligibility_conflict)) {
    stop(
      "cannabis_current and cannabis_frequency are not complementary: ",
      sum(eligibility_conflict), " conflicting rows.",
      call. = FALSE
    )
  }

  audit_cudit_zero_rule <- function(data, dataset_label, wave_column = NULL) {
    audit_one <- function(subset, wave_label) {
      subset_frequency <- suppressWarnings(as.numeric(subset[["cannabis_frequency"]]))
      subset_current <- suppressWarnings(as.numeric(subset[[current_use_var]]))
      confirmed_nonuser <- !is.na(subset_frequency) & !is.na(subset_current) &
        subset_current == 0 & subset_frequency == 0
      confirmed_user <- !is.na(subset_frequency) & !is.na(subset_current) &
        subset_current == 1 & subset_frequency > 0
      data.table::rbindlist(lapply(c("cudit_freq", "cudit_sev"), function(variable) {
        invalid_nonuser <- confirmed_nonuser &
          (is.na(subset[[variable]]) | subset[[variable]] != 0)
        if (any(invalid_nonuser)) {
          stop(
            dataset_label, " ", wave_label, ": ", variable,
            " must be zero for every confirmed non-user; invalid rows: ",
            sum(invalid_nonuser), ". Rebuild with the corrected derived_data_v9.",
            call. = FALSE
          )
        }
        data.table::data.table(
          dataset = dataset_label,
          wave = wave_label,
          variable = variable,
          confirmed_nonusers = sum(confirmed_nonuser),
          confirmed_nonusers_zero = sum(subset[[variable]][confirmed_nonuser] == 0),
          confirmed_users = sum(confirmed_user),
          confirmed_user_missing = sum(is.na(subset[[variable]][confirmed_user]))
        )
      }))
    }
    if (is.null(wave_column)) return(audit_one(data, "baseline"))
    data.table::rbindlist(lapply(unique(data[[wave_column]]), function(wave_label) {
      audit_one(data[get(wave_column) == wave_label], as.character(wave_label))
    }))
  }

  cudit_zero_audit <- data.table::rbindlist(list(
    audit_cudit_zero_rule(baseline, "baseline_genomic"),
    audit_cudit_zero_rule(baseline_longitudinal, "baseline_longitudinal"),
    audit_cudit_zero_rule(longitudinal, "longitudinal", wave_var)
  ))
  write_csv(
    cudit_zero_audit,
    file.path(output_dir, "qc", "cudit_structural_zero_audit.csv")
  )

  zero_coded_variables <- phenotype_manifest[
    primary_population == "all" & zero_coded_sensitivity == TRUE,
    variable
  ]
  audit_zero_coded_outcomes <- function(
      data,
      dataset_label,
      variables,
      wave_label = "baseline") {
    confirmed_nonuser <- !is.na(data[[current_use_var]]) &
      !is.na(data[["cannabis_frequency"]]) &
      data[[current_use_var]] == 0 & data[["cannabis_frequency"]] == 0
    variables <- intersect(variables, names(data))
    data.table::rbindlist(lapply(variables, function(variable) {
      invalid <- confirmed_nonuser &
        (is.na(data[[variable]]) | data[[variable]] != 0)
      if (any(invalid)) {
        stop(
          dataset_label, " ", wave_label, ": ", variable,
          " must be zero for every confirmed non-user; invalid rows: ",
          sum(invalid), ".",
          call. = FALSE
        )
      }
      data.table::data.table(
        dataset = dataset_label,
        wave = wave_label,
        variable = variable,
        confirmed_nonusers = sum(confirmed_nonuser),
        confirmed_nonusers_zero = sum(data[[variable]][confirmed_nonuser] == 0),
        confirmed_users = sum(
          !is.na(data[[current_use_var]]) & data[[current_use_var]] == 1
        ),
        confirmed_user_missing = sum(
          !is.na(data[[current_use_var]]) & data[[current_use_var]] == 1 &
            is.na(data[[variable]])
        )
      )
    }))
  }
  zero_coded_audit <- data.table::rbindlist(list(
    audit_zero_coded_outcomes(
      baseline, "baseline_genomic", zero_coded_variables
    ),
    audit_zero_coded_outcomes(
      baseline_longitudinal,
      "baseline_longitudinal",
      zero_coded_variables
    ),
    data.table::rbindlist(lapply(unique(longitudinal[[wave_var]]), function(wave_label) {
      audit_zero_coded_outcomes(
        longitudinal[get(wave_var) == wave_label],
        "longitudinal",
        zero_coded_variables,
        as.character(wave_label)
      )
    }))
  ), fill = TRUE)
  write_csv(
    zero_coded_audit,
    file.path(output_dir, "qc", "structural_zero_outcome_audit.csv")
  )

  cud_status_zero_audit <- data.table::rbindlist(lapply(
    unique(longitudinal[[wave_var]]),
    function(wave_label) {
      subset <- longitudinal[get(wave_var) == wave_label]
      confirmed_nonuser <- !is.na(subset[[current_use_var]]) &
        !is.na(subset[["cannabis_frequency"]]) &
        subset[[current_use_var]] == 0 & subset[["cannabis_frequency"]] == 0
      invalid <- confirmed_nonuser &
        (
          is.na(subset[[cud_var]]) | subset[[cud_var]] != 0 |
            is.na(subset[[cud_symptom_var]]) |
            subset[[cud_symptom_var]] != 0
        )
      if (any(invalid)) {
        stop(
          "Longitudinal ", wave_label,
          ": confirmed same-wave non-users must have CUD status zero; invalid rows: ",
          sum(invalid), ".",
          call. = FALSE
        )
      }
      data.table::data.table(
        wave = as.character(wave_label),
        rows = nrow(subset),
        cud_observed = sum(!is.na(subset[[cud_var]])),
        cud_cases = sum(subset[[cud_var]] == 1, na.rm = TRUE),
        confirmed_nonusers = sum(confirmed_nonuser),
        confirmed_nonusers_cud_zero = sum(subset[[cud_var]][confirmed_nonuser] == 0),
        confirmed_nonusers_symptom_zero = sum(
          subset[[cud_symptom_var]][confirmed_nonuser] == 0
        ),
        complete_symptom_assessments = sum(
          subset[[cud_symptom_observed_var]] == 11,
          na.rm = TRUE
        ),
        partial_symptom_assessments = sum(
          subset[[cud_symptom_observed_var]] > 0 &
            subset[[cud_symptom_observed_var]] < 11,
          na.rm = TRUE
        ),
        unresolved_cud_missing = sum(
          is.na(subset[[current_use_var]]) & is.na(subset[[cud_var]])
        )
      )
    }
  ))
  write_csv(
    cud_status_zero_audit,
    file.path(output_dir, "qc", "cud_status_structural_zero_audit.csv")
  )

  range_audit <- phenotype_manifest[, {
    x <- baseline[[variable]]
    observed_min <- suppressWarnings(min(as.numeric(x), na.rm = TRUE))
    observed_max <- suppressWarnings(max(as.numeric(x), na.rm = TRUE))
    if (!is.finite(observed_min)) observed_min <- NA_real_
    if (!is.finite(observed_max)) observed_max <- NA_real_
    expected_min <- suppressWarnings(as.numeric(valid_min))
    expected_max <- suppressWarnings(as.numeric(valid_max))
    data.table::data.table(
      n_nonmissing = sum(!is.na(x)),
      observed_min = observed_min,
      observed_max = observed_max,
      expected_min = expected_min,
      expected_max = expected_max,
      below_expected = ifelse(is.na(expected_min), FALSE, observed_min < expected_min),
      above_expected = ifelse(is.na(expected_max), FALSE, observed_max > expected_max)
    )
  }, by = .(variable, label, family)]
  write_csv(range_audit, file.path(output_dir, "qc", "phenotype_range_audit.csv"))
  if (any(range_audit$below_expected | range_audit$above_expected, na.rm = TRUE)) {
    stop(
      "Phenotype range audit failed. Review qc/phenotype_range_audit.csv before analysis.",
      call. = FALSE
    )
  }

  pgs_qc <- pgs_manifest[, {
    x <- as.numeric(baseline[[variable]])
    data.table::data.table(
      n = length(x),
      n_missing = sum(is.na(x)),
      missing_percent = 100 * mean(is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }, by = .(variable, label, domain)]
  write_csv(pgs_qc, file.path(output_dir, "qc", "pgs_distribution_qc.csv"))

  list(
    baseline = baseline,
    baseline_longitudinal = baseline_longitudinal,
    longitudinal = longitudinal,
    pgs_manifest = pgs_manifest,
    phenotype_manifest = phenotype_manifest,
    coding_audit = coding,
    range_audit = range_audit,
    pgs_qc = pgs_qc
  )
}
