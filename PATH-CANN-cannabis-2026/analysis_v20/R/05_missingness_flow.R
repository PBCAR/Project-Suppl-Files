summarize_v20_longitudinal_missingness <- function(
    longitudinal,
    wave_var,
    variables,
    waves) {
  data.table::rbindlist(lapply(waves, function(wave) {
    wave_data <- longitudinal[get(wave_var) == wave]
    data.table::rbindlist(lapply(variables, function(variable) {
      data.table::data.table(
        wave = wave,
        variable = variable,
        rows = nrow(wave_data),
        observed = sum(!is.na(wave_data[[variable]])),
        missing = sum(is.na(wave_data[[variable]])),
        missing_percent = 100 * mean(is.na(wave_data[[variable]]))
      )
    }))
  }))
}

make_v20_participant_flow <- function(
    baseline,
    baseline_longitudinal,
    longitudinal,
    id_var,
    wave_var,
    waves,
    cud_var,
    cud_symptom_observed_var,
    min_complete_symptoms) {
  wave_rows <- data.table::rbindlist(lapply(waves, function(wave) {
    wave_data <- longitudinal[get(wave_var) == wave]
    data.table::data.table(
      stage = c(
        paste0(wave, " CUD status observed"),
        paste0(wave, " complete 11-item symptom assessment")
      ),
      participants = c(
        sum(!is.na(wave_data[[cud_var]])),
        sum(
          wave_data[[cud_symptom_observed_var]] >= min_complete_symptoms,
          na.rm = TRUE
        )
      ),
      cud_positive = c(
        sum(wave_data[[cud_var]] == 1, na.rm = TRUE),
        sum(
          wave_data[[cud_var]] == 1 &
            wave_data[[cud_symptom_observed_var]] >= min_complete_symptoms,
          na.rm = TRUE
        )
      )
    )
  }))
  data.table::rbindlist(list(
    data.table::data.table(
      stage = c("Baseline full phenotype cohort", "Baseline genomic cohort"),
      participants = c(
        data.table::uniqueN(baseline_longitudinal[[id_var]]),
        data.table::uniqueN(baseline[[id_var]])
      ),
      cud_positive = NA_integer_
    ),
    wave_rows
  ), fill = TRUE)
}

run_v20_missingness_flow <- function(
    baseline,
    baseline_longitudinal,
    longitudinal,
    id_var,
    wave_var,
    waves,
    variables,
    cud_var,
    cud_symptom_observed_var,
    min_complete_symptoms,
    output_dir) {
  missingness <- summarize_v20_longitudinal_missingness(
    longitudinal, wave_var, variables, waves
  )
  flow <- make_v20_participant_flow(
    baseline, baseline_longitudinal, longitudinal, id_var, wave_var, waves,
    cud_var, cud_symptom_observed_var, min_complete_symptoms
  )
  write_csv(
    missingness,
    file.path(output_dir, "missingness", "longitudinal_variable_missingness.csv")
  )
  write_csv(
    flow,
    file.path(output_dir, "missingness", "participant_flow.csv")
  )
  list(missingness = missingness, flow = flow)
}
