mean_ci <- function(x) {
  observed <- as.numeric(x[!is.na(x)])
  n <- length(observed)
  if (!n) {
    return(data.table::data.table(
      n = 0L, estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_
    ))
  }
  estimate <- mean(observed)
  if (n < 2L) {
    return(data.table::data.table(
      n = n, estimate = estimate, conf_low = NA_real_, conf_high = NA_real_
    ))
  }
  standard_error <- stats::sd(observed) / sqrt(n)
  critical <- stats::qt(0.975, df = n - 1L)
  data.table::data.table(
    n = n,
    estimate = estimate,
    conf_low = estimate - critical * standard_error,
    conf_high = estimate + critical * standard_error
  )
}

proportion_ci <- function(x) {
  observed <- as.numeric(x[!is.na(x)])
  n <- length(observed)
  if (!n) {
    return(data.table::data.table(
      n = 0L, estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_
    ))
  }
  estimate <- mean(observed)
  interval <- suppressWarnings(stats::prop.test(
    sum(observed == 1), n, correct = FALSE
  )$conf.int)
  data.table::data.table(
    n = n,
    estimate = estimate,
    conf_low = interval[[1]],
    conf_high = interval[[2]]
  )
}

summarize_longitudinal_measure <- function(
    data,
    wave,
    measure,
    label,
    population,
    scale_max,
    summary_type = "mean") {
  summary <- if (summary_type == "proportion") {
    proportion_ci(data[[measure]])
  } else {
    mean_ci(data[[measure]])
  }
  cbind(
    data.table::data.table(
      wave = wave,
      wave_number = wave_number(wave),
      measure = measure,
      measure_label = label,
      population = population,
      scale_max = scale_max,
      summary_type = summary_type
    ),
    summary
  )
}
run_longitudinal_descriptives <- function(
    baseline,
    longitudinal,
    display_waves,
    id_var,
    wave_var,
    current_use_var,
    cud_var,
    cud_symptom_var,
    cud_symptom_observed_var,
    min_complete_symptoms,
    output_dir) {
  
  baseline_required <- c(
    id_var,
    current_use_var,
    "cudit_freq",
    "cudit_sev"
  )
  
  longitudinal_required <- c(
    id_var,
    wave_var,
    current_use_var,
    "cudit_freq",
    "cudit_sev",
    cud_var,
    cud_symptom_var,
    cud_symptom_observed_var
  )
  
  check_required_columns(
    baseline,
    baseline_required,
    "Baseline descriptive data"
  )
  
  check_required_columns(
    longitudinal,
    longitudinal_required,
    "Longitudinal descriptive data"
  )
  
  # Extract the longitudinal wave variable explicitly. This avoids ambiguity
  # between the name of the data column and the value used by the loop.
  longitudinal_wave_values <- as.character(
    longitudinal[[wave_var]]
  )
  
  requested_followup_waves <- setdiff(
    as.character(display_waves),
    "T1"
  )
  
  available_followup_waves <- unique(
    longitudinal_wave_values[
      !is.na(longitudinal_wave_values)
    ]
  )
  
  missing_followup_waves <- setdiff(
    requested_followup_waves,
    available_followup_waves
  )
  
  if (length(missing_followup_waves)) {
    stop(
      paste0(
        "The longitudinal data do not contain requested wave(s): ",
        paste(missing_followup_waves, collapse = ", "),
        ". Available waves are: ",
        paste(sort(available_followup_waves), collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  summaries <- list()
  index <- 1L
  
  for (wave_value in as.character(display_waves)) {
    
    if (wave_value == "T1") {
      wave_data <- data.table::copy(
        baseline[, ..baseline_required]
      )
      
      wave_data[, `:=`(
        cud_case_plot = NA_real_,
        cud_symptom_plot = NA_real_,
        symptom_complete_plot = FALSE
      )]
      
    } else {
      # Standard-evaluation subsetting: longitudinal_wave_values is an
      # ordinary vector and wave_value is the current loop value.
      wave_rows <- which(
        longitudinal_wave_values == wave_value
      )
      
      if (!length(wave_rows)) {
        stop(
          paste0(
            "No longitudinal observations were found for ",
            wave_var,
            " = ",
            wave_value,
            "."
          ),
          call. = FALSE
        )
      }
      
      wave_data <- data.table::copy(
        longitudinal[wave_rows]
      )
      
      # Confirm that the subset contains exactly the intended wave.
      observed_wave_values <- unique(
        as.character(wave_data[[wave_var]])
      )
      
      observed_wave_values <- observed_wave_values[
        !is.na(observed_wave_values)
      ]
      
      if (
        length(observed_wave_values) != 1L ||
        observed_wave_values[[1]] != wave_value
      ) {
        stop(
          paste0(
            "Wave subsetting failed for ",
            wave_value,
            ". Observed wave values were: ",
            paste(observed_wave_values, collapse = ", "),
            "."
          ),
          call. = FALSE
        )
      }
      
      wave_data[, `:=`(
        cud_case_plot = as.numeric(
          get(cud_var)
        ),
        cud_symptom_plot = as.numeric(
          get(cud_symptom_var)
        ),
        symptom_complete_plot =
          as.numeric(
            get(cud_symptom_observed_var)
          ) >= min_complete_symptoms
      )]
    }
    
    summaries[[index]] <- summarize_longitudinal_measure(
      wave_data,
      wave_value,
      current_use_var,
      "Current cannabis use",
      "all",
      1,
      "proportion"
    )
    index <- index + 1L
    
    if (wave_value != "T1") {
      summaries[[index]] <- summarize_longitudinal_measure(
        wave_data,
        wave_value,
        "cud_case_plot",
        "DSM-5 CUD",
        "all",
        1,
        "proportion"
      )
      index <- index + 1L
    }
    
    for (spec in list(
      list(
        variable = "cudit_freq",
        label = "CUDIT-C",
        maximum = 8
      ),
      list(
        variable = "cudit_sev",
        label = "CUDIT-P",
        maximum = 24
      )
    )) {
      summaries[[index]] <- summarize_longitudinal_measure(
        wave_data,
        wave_value,
        spec$variable,
        spec$label,
        "all",
        spec$maximum,
        "mean"
      )
      index <- index + 1L
    }
    
    if (wave_value != "T1") {
      symptom_all <- wave_data[
        (
          !is.na(get(current_use_var)) &
            get(current_use_var) == 0 &
            cud_symptom_plot == 0
        ) |
          symptom_complete_plot == TRUE
      ]
      
      summaries[[index]] <- summarize_longitudinal_measure(
        symptom_all,
        wave_value,
        "cud_symptom_plot",
        "DSM-5 CUD symptoms",
        "all_eligible_symptoms",
        11,
        "mean"
      )
      index <- index + 1L
    }
    
    current_users <- wave_data[
      !is.na(get(current_use_var)) &
        get(current_use_var) == 1
    ]
    
    for (spec in list(
      list(
        variable = "cudit_freq",
        label = "CUDIT-C",
        maximum = 8
      ),
      list(
        variable = "cudit_sev",
        label = "CUDIT-P",
        maximum = 24
      )
    )) {
      summaries[[index]] <- summarize_longitudinal_measure(
        current_users,
        wave_value,
        spec$variable,
        spec$label,
        "current_users",
        spec$maximum,
        "mean"
      )
      index <- index + 1L
    }
    
    if (wave_value != "T1") {
      symptom_users <- current_users[
        symptom_complete_plot == TRUE
      ]
      
      summaries[[index]] <- summarize_longitudinal_measure(
        symptom_users,
        wave_value,
        "cud_symptom_plot",
        "DSM-5 CUD symptoms",
        "current_users_complete_symptoms",
        11,
        "mean"
      )
      index <- index + 1L
    }
  }
  
  result <- data.table::rbindlist(
    summaries,
    fill = TRUE
  )
  
  result[
    ,
    estimate_percent_of_scale :=
      100 * estimate / scale_max
  ]
  
  result[
    ,
    conf_low_percent_of_scale :=
      100 * conf_low / scale_max
  ]
  
  result[
    ,
    conf_high_percent_of_scale :=
      100 * conf_high / scale_max
  ]
  
  result[
    summary_type == "proportion",
    `:=`(
      estimate_percent_of_scale =
        100 * estimate,
      conf_low_percent_of_scale =
        100 * conf_low,
      conf_high_percent_of_scale =
        100 * conf_high
    )
  ]
  
  data.table::setorder(
    result,
    wave_number,
    population,
    measure
  )
  
  write_csv(
    result,
    file.path(
      output_dir,
      "descriptives",
      "longitudinal_cannabis_profile.csv"
    )
  )
  
  result
}