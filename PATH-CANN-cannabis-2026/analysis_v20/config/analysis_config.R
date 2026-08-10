analysis_root <- normalizePath(
  Sys.getenv("PATH_CANN_ANALYSIS_ROOT", unset = getwd()),
  mustWork = FALSE
)

if (!file.exists(file.path(analysis_root, "config", "phenotype_manifest.csv"))) {
  candidate <- normalizePath(file.path(getwd(), "analysis_v20"), mustWork = FALSE)
  if (file.exists(file.path(candidate, "config", "phenotype_manifest.csv"))) {
    analysis_root <- candidate
  }
}

storage_root <- normalizePath(
  Sys.getenv(
    "PATH_CANN_STORAGE_ROOT",
    unset = file.path(analysis_root, "..", "..")
  ),
  mustWork = FALSE
)

input_dir <- "/Users/weideng/Library/CloudStorage/OneDrive-McMasterUniversity/PBCAR/PATH/cannabis-2025/"

baseline_file <- Sys.getenv(
  "PATH_CANN_BASELINE_FILE",
  unset = file.path(input_dir, "baseline_analysis.rds")
)

baseline_longitudinal_file <- Sys.getenv(
  "PATH_CANN_BASELINE_LONGITUDINAL_FILE",
  unset = file.path(input_dir, "baseline_longitudinal_predictors.rds")
)

longitudinal_file <- Sys.getenv(
  "PATH_CANN_LONGITUDINAL_FILE",
  unset = file.path(input_dir, "longitudinal_analysis.rds")
)

run_id <- Sys.getenv(
  "PATH_CANN_RUN_ID",
  unset = paste0(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "_v20")
)

output_base <- normalizePath(
  Sys.getenv(
    "PATH_CANN_OUTPUT_BASE",
    unset = file.path(
      storage_root, "PATH_CANN_results", "analysis_v20", "outputs"
    )
  ),
  mustWork = FALSE
)
output_dir <- file.path(output_base, run_id)
manifest_dir <- file.path(analysis_root, "config")
pgs_manifest_file <- Sys.getenv(
  "PATH_CANN_PGS_MANIFEST",
  unset = file.path(manifest_dir, "pgs_manifest.csv")
)
phenotype_manifest_file <- Sys.getenv(
  "PATH_CANN_PHENOTYPE_MANIFEST",
  unset = file.path(manifest_dir, "phenotype_manifest.csv")
)
table_characteristics_manifest_file <- Sys.getenv(
  "PATH_CANN_TABLE_CHARACTERISTICS_MANIFEST",
  unset = file.path(manifest_dir, "table_characteristics_manifest.csv")
)

id_var <- "id"
baseline_age_var <- "age_t1"
wave_age_var <- "age_wave"
sex_var <- "sex"
wave_var <- "wave"
cud_var <- "cud_case"
cud_symptom_var <- "cud_symptom_count"
cud_symptom_observed_var <- "cud_symptom_n_observed"
current_use_var <- "cannabis_current"
pc_vars <- paste0("PC", 1:10)

baseline_frequency_var <- "cudit_freq"
baseline_coping_var <- "mmm_coping"
baseline_severity_var <- "cudit_sev"
baseline_perceived_harm_var <- "nsduh3"

concurrent_frequency_var <- "cudit_freq"
concurrent_coping_var <- "mmm_coping"
concurrent_severity_var <- "cudit_sev"
concurrent_perceived_harm_var <- "nsduh3"

prediction_waves <- paste0("T", 12:17)
display_waves <- c("T1", prediction_waves)

# The wave-specific CUD analysis selects PGSs from the baseline association
# map: externalizing and CUD are retained as a priori anchors, and another PGS
# is retained when it has at least one global-FDR phenotype outside the
# externalizing-associated phenotype set. The selection audit is written with
# every run. Wave-specific estimates and omnibus interaction tests are
# marginal; no repeated-measures mixed model is fitted.
cud_followup_anchor_pgs <- c("pgs_externalizing", "pgs_cud")
complete_cud_symptoms_required <- 11L

seed <- as.integer(Sys.getenv("PATH_CANN_SEED", unset = "20260721"))

global_fdr_alpha <- 0.05
