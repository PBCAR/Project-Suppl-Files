# PATH-CANN analysis v20 production entry point.
# Open this file in RStudio and click Source, or run it with Rscript.

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_root <- if (length(file_argument)) {
  dirname(normalizePath(sub("^--file=", "", file_argument[[1]]), mustWork = FALSE))
} else if (
  requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
    nzchar(rstudioapi::getSourceEditorContext()$path)
) {
  dirname(normalizePath(rstudioapi::getSourceEditorContext()$path, mustWork = FALSE))
} else {
  normalizePath(getwd(), mustWork = FALSE)
}

if (!file.exists(file.path(script_root, "config", "analysis_config.R"))) {
  stop("Open analysis_v20/run_all.R in RStudio and click Source.", call. = FALSE)
}

Sys.setenv(PATH_CANN_ANALYSIS_ROOT = script_root)
source(file.path(script_root, "config", "analysis_config.R"))

production_modules <- c(
  "00_utils.R",
  "01_validate_inputs.R",
  "02_baseline_pgs_screen.R",
  "04_longitudinal_descriptives.R",
  "05_cud_pgs_selection.R",
  "05_cud_pgs_marginal.R",
  "05_missingness_flow.R",
  "05_figures.R",
  "06_participant_tables.R"
)
production_paths <- file.path(analysis_root, "R", production_modules)
missing_modules <- production_modules[!file.exists(production_paths)]
if (length(missing_modules)) {
  stop(
    "Missing production module(s): ", paste(missing_modules, collapse = ", "),
    call. = FALSE
  )
}
for (file in production_paths) source(file)

check_packages()
config_snapshot <- c(
  paste0("Baseline genomic file: ", baseline_file),
  paste0("Full phenotype baseline file: ", baseline_longitudinal_file),
  paste0("Longitudinal file: ", longitudinal_file),
  "Baseline phenotype coding: medical/recreational no-opinion responses are missing",
  "MPT coding: log1p of each raw non-negative parameter, recomputed at run time",
  "All-participant structural-zero outcomes: CUDIT-C, CUDIT-P, all MMM motives, and expenditure",
  "Current-user outcomes: MPT and cannabis-effect measures",
  paste0("Repeated CUD waves: ", paste(prediction_waves, collapse = ", ")),
  paste0("CUD follow-up anchor PGSs: ", paste(cud_followup_anchor_pgs, collapse = ", ")),
  paste0(
    "CUD follow-up selection: both anchors plus any other PGS with at least ",
    "one baseline global-FDR phenotype outside the externalizing hit set"
  ),
  "CUD follow-up model: marginal logistic model with wave-saturated adjustment",
  "Wave-specific CUD estimates: marginal odds ratios with participant-clustered HC1 covariance",
  "Temporal heterogeneity: omnibus PGS-by-wave Wald test using the clustered covariance",
  "No random-intercept mixed, phenotype-adjustment, mediation, attenuation, AUC, or clinical prediction model",
  "Figure S3 cohort: baseline genetic analysis participants only",
  "Participant characteristics and all-phenotype distributions: supplementary tables",
  "Table S2 sex comparisons: designated-population tests with BH q-values across 30 phenotypes",
  paste0("Global FDR alpha: ", global_fdr_alpha),
  "Baseline FDR family: 360 designated-population PGS-phenotype tests",
  "Repeated-CUD FDR families: selected-PGS-by-wave estimates and selected-PGS omnibus interactions"
)

initialize_output(output_dir, run_id, config_snapshot)

inputs <- load_and_validate_inputs(
  baseline_file, baseline_longitudinal_file, longitudinal_file,
  pgs_manifest_file, phenotype_manifest_file,
  id_var, baseline_age_var, wave_age_var, sex_var, wave_var, cud_var,
  cud_symptom_var, cud_symptom_observed_var,
  current_use_var, pc_vars, output_dir
)

covariates <- c(baseline_age_var, sex_var, pc_vars)
if ("genotyping_batch" %in% names(inputs$baseline)) {
  covariates <- c(covariates, "genotyping_batch")
}

message("Running the locked 12-PGS by 30-phenotype baseline screen...")
baseline_results <- run_baseline_pgs_screen(
  inputs$baseline, inputs$pgs_manifest, inputs$phenotype_manifest,
  covariates, current_use_var, output_dir, global_fdr_alpha
)

message("Creating supplementary participant-characteristics and phenotype tables...")
participant_tables <- run_baseline_participant_tables(
  inputs$baseline, inputs$phenotype_manifest,
  table_characteristics_manifest_file,
  current_use_var, sex_var, output_dir
)

message("Selecting PGSs for repeated CUD follow-up from the baseline association map...")
cud_pgs_selection <- select_cud_followup_pgs(
  baseline_results$primary, inputs$pgs_manifest, output_dir,
  externalizing_pgs = cud_followup_anchor_pgs[[1]],
  cud_pgs = cud_followup_anchor_pgs[[2]]
)

message("Running selected-PGS marginal CUD models and omnibus PGS-by-wave Wald tests...")
cud_pgs_results <- run_selected_cud_pgs_marginal(
  inputs$baseline, inputs$longitudinal,
  cud_pgs_selection$audit, prediction_waves,
  id_var, sex_var, wave_var, wave_age_var, cud_var,
  pc_vars, output_dir, global_fdr_alpha
)
cud_pgs_results$selection <- cud_pgs_selection$audit

message("Summarizing Figure S3 in the baseline genetic analysis cohort...")
genetic_ids <- unique(inputs$baseline[[id_var]])
genetic_longitudinal <- data.table::copy(
  inputs$longitudinal[get(id_var) %in% genetic_ids]
)
if (any(!genetic_longitudinal[[id_var]] %in% genetic_ids)) {
  stop("Figure S3 contains a follow-up ID outside the baseline genetic cohort.", call. = FALSE)
}
longitudinal_profile <- run_longitudinal_descriptives(
  inputs$baseline, genetic_longitudinal, display_waves,
  id_var, wave_var, current_use_var, cud_var, cud_symptom_var,
  cud_symptom_observed_var, complete_cud_symptoms_required, output_dir
)
figure_s3_followup_audit <- genetic_longitudinal[, .(
  genetic_rows = .N,
  genetic_participants = data.table::uniqueN(get(id_var)),
  cud_observed = sum(!is.na(get(cud_var))),
  cud_cases = sum(get(cud_var) == 1, na.rm = TRUE),
  cud_prevalence = mean(get(cud_var), na.rm = TRUE)
), by = c(wave_var)]
figure_s3_full_audit <- inputs$longitudinal[, .(
  full_rows = .N,
  full_participants = data.table::uniqueN(get(id_var))
), by = c(wave_var)]
figure_s3_followup_audit <- merge(
  figure_s3_followup_audit,
  figure_s3_full_audit,
  by = wave_var,
  all.x = TRUE
)
cud_profile_audit <- longitudinal_profile[
  measure == "cud_case_plot" & population == "all",
  .(
    wave,
    profile_cud_denominator = n,
    profile_cud_prevalence = estimate
  )
]
figure_s3_followup_audit <- merge(
  figure_s3_followup_audit,
  cud_profile_audit,
  by.x = wave_var,
  by.y = "wave",
  all.x = TRUE
)
figure_s3_followup_audit[, `:=`(
  excluded_non_genetic_participants = full_participants - genetic_participants,
  profile_matches_genetic_subset =
    cud_observed == profile_cud_denominator &
      abs(cud_prevalence - profile_cud_prevalence) < 1e-12
)]
if (any(!figure_s3_followup_audit$profile_matches_genetic_subset)) {
  stop(
    "Figure S3 CUD prevalence does not match the genetic-cohort wave audit.",
    call. = FALSE
  )
}
figure_s3_cohort_audit <- data.table::rbindlist(list(
  data.table::data.table(
    wave = "T1",
    genetic_rows = nrow(inputs$baseline),
    genetic_participants = data.table::uniqueN(inputs$baseline[[id_var]]),
    cud_observed = NA_integer_,
    cud_cases = NA_integer_,
    cud_prevalence = NA_real_,
    full_rows = nrow(inputs$baseline_longitudinal),
    full_participants = data.table::uniqueN(inputs$baseline_longitudinal[[id_var]]),
    profile_cud_denominator = NA_integer_,
    profile_cud_prevalence = NA_real_,
    excluded_non_genetic_participants =
      data.table::uniqueN(inputs$baseline_longitudinal[[id_var]]) -
        data.table::uniqueN(inputs$baseline[[id_var]]),
    profile_matches_genetic_subset = TRUE
  ),
  figure_s3_followup_audit
), use.names = TRUE, fill = TRUE)
data.table::setnames(figure_s3_cohort_audit, wave_var, "wave", skip_absent = TRUE)
write_csv(
  figure_s3_cohort_audit,
  file.path(output_dir, "qc", "figure_s3_genetic_cohort_audit.csv")
)

message("Writing missingness and participant-flow outputs...")
missingness_results <- run_v20_missingness_flow(
  inputs$baseline, inputs$baseline_longitudinal, inputs$longitudinal,
  id_var, wave_var, prediction_waves,
  c(
    current_use_var, "cudit_freq", "cudit_sev", "mmm_coping",
    cud_var, cud_symptom_var, cud_symptom_observed_var
  ),
  cud_var, cud_symptom_observed_var, complete_cud_symptoms_required,
  output_dir
)

message("Creating v20 figures and evidence table...")
figure_results <- run_figures(
  baseline_results, cud_pgs_results, longitudinal_profile, output_dir
)

summary_lines <- c(
  "PATH-CANN analysis revision v20 completed successfully.",
  "CUD follow-up estimand: selected-PGS marginal association with DSM-5 CUD status at T12-T17.",
  "All-participant structural-zero outcomes: CUDIT-C, CUDIT-P, MMM motives, and expenditure.",
  "Current-user outcomes: MPT and cannabis-effect measures.",
  paste0("Baseline primary tests: ", nrow(baseline_results$primary), "."),
  paste0("Baseline global-FDR hits: ", nrow(baseline_results$significant), "."),
  paste0(
    "No-opinion responses recoded missing: ",
    sum(inputs$coding_audit$acceptability_audit$no_opinion_recode_to_missing), "."
  ),
  paste0(
    "PGSs selected for CUD follow-up: ",
    paste(cud_pgs_selection$manifest$label, collapse = ", "), "."
  ),
  paste0(
    "Wave-specific marginal CUD FDR hits: ",
    sum(cud_pgs_results$wave_effects$passes_wave_fdr, na.rm = TRUE), "."
  ),
  paste0(
    "PGS-by-wave omnibus Wald FDR hits: ",
    sum(cud_pgs_results$interactions$passes_interaction_fdr, na.rm = TRUE), "."
  ),
  paste0(
    "Figure S3 baseline genetic cohort N: ",
    data.table::uniqueN(inputs$baseline[[id_var]]), "."
  ),
  "Figure S3 follow-up rows are restricted to IDs in that baseline genetic cohort.",
  "No random-intercept mixed, AUC, screening, clinical prediction, phenotype-adjustment, mediation, or attenuation model was run.",
  "Differences in wave-specific significance do not by themselves demonstrate temporal heterogeneity.",
  paste0("Output directory: ", output_dir)
)
writeLines(
  summary_lines,
  file.path(output_dir, "analysis_summary.txt"),
  useBytes = TRUE
)
cat(paste(summary_lines, collapse = "\n"), "\n")
