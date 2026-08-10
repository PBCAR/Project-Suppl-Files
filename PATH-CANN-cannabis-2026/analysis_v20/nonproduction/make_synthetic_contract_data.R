args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else file.path(tempdir(), "path_cann_v20_synthetic")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260731)
n <- 450L
ids <- sprintf("S%04d", seq_len(n))
pc_vars <- paste0("PC", 1:10)
pgs_manifest <- data.table::fread(file.path(
  Sys.getenv("PATH_CANN_ANALYSIS_ROOT"), "config", "pgs_manifest.csv"
))
phenotype_manifest <- data.table::fread(file.path(
  Sys.getenv("PATH_CANN_ANALYSIS_ROOT"), "config", "phenotype_manifest.csv"
))

latent_externalizing <- rnorm(n)
latent_somatic <- rnorm(n)
latent_adhd <- 0.45 * latent_externalizing + rnorm(n, sd = 0.9)
random_intercept <- rnorm(n, sd = 0.65)

baseline <- data.table::data.table(
  id = ids,
  age_t1 = runif(n, 25, 65),
  sex = sample(c("Female", "Male"), n, replace = TRUE),
  ethnicity = rbinom(n, 1, 0.22),
  financial2 = sample(1:9, n, replace = TRUE),
  genotyping_batch = sample(c("A", "B", "C"), n, replace = TRUE)
)
for (pc in pc_vars) baseline[, (pc) := rnorm(n)]
for (pgs in pgs_manifest$variable) baseline[, (pgs) := rnorm(n)]
baseline[, pgs_externalizing := latent_externalizing]
baseline[, pgs_adhd := latent_adhd]
baseline[, pgs_somatic_factor := latent_somatic]

baseline[, cannabis_ever := rbinom(
  n, 1, plogis(0.7 + 0.35 * latent_externalizing)
)]
baseline[, cannabis_current := rbinom(
  n, 1, plogis(-0.05 + 0.45 * latent_externalizing)
)]
baseline[cannabis_current == 1, cannabis_ever := 1]
baseline[, cannabis_frequency := ifelse(
  cannabis_current == 1, sample(1:7, .N, replace = TRUE), 0
)]
baseline[, cannabis_firstuse := ifelse(
  cannabis_ever == 1, pmax(12, round(rnorm(.N, 19, 3))), NA_real_
)]
baseline[, cudit_freq := ifelse(
  cannabis_current == 1,
  pmin(8, rpois(.N, lambda = exp(0.8 + 0.20 * latent_externalizing))),
  0
)]
baseline[, cudit_sev := ifelse(
  cannabis_current == 1,
  pmin(24, rpois(.N, lambda = exp(0.55 + 0.18 * cudit_freq +
    0.16 * latent_adhd))),
  0
)]

zero_continuous <- c(
  "cannabis_expenditure", "mmm_enhancement", "mmm_conformity",
  "mmm_expansion", "mmm_coping", "mmm_social"
)
for (variable in zero_continuous) {
  baseline[, (variable) := ifelse(
    cannabis_current == 1,
    if (variable == "cannabis_expenditure") {
      pmax(0, rgamma(.N, shape = 2, scale = 25))
    } else {
      pmin(20, pmax(0, rnorm(.N, 7 + 0.8 * latent_externalizing, 3)))
    },
    0
  )]
}
baseline[, mmm_coping := ifelse(
  cannabis_current == 1,
  pmin(20, pmax(0, rnorm(.N, 7 + latent_adhd, 3))),
  0
)]

for (variable in c("Breakpoint", "Intensity", "Omax", "Pmax", "Alpha")) {
  baseline[, (variable) := ifelse(
    cannabis_current == 1,
    pmax(0, rgamma(.N, shape = 1.5, scale = 2)),
    NA_real_
  )]
  baseline[, (paste0("log", variable)) := NA_real_]
}
for (variable in c("cann_med_acceptable", "cann_rec_acceptable", "cannabis_acceptable")) {
  baseline[, (variable) := sample(1:5, .N, replace = TRUE)]
}
baseline[, nsduh3 := sample(1:4, .N, replace = TRUE)]
baseline[, cannabis_risk := pmin(12, rpois(.N, 4))]
baseline[, cannabis_benefits := pmin(11, rpois(.N, 4))]
for (variable in paste0("cannabis_impact", 1:7)) {
  baseline[, (variable) := ifelse(
    cannabis_current == 1, sample(1:5, .N, replace = TRUE), NA_integer_
  )]
}
baseline[, mae := ifelse(
  cannabis_current == 1, pmin(70, pmax(14, round(rnorm(.N, 35, 8)))), NA_real_
)]

missing_phenotypes <- setdiff(phenotype_manifest$variable, names(baseline))
if (length(missing_phenotypes)) {
  stop("Synthetic generator omitted phenotype(s): ", paste(missing_phenotypes, collapse = ", "))
}

baseline_longitudinal <- data.table::copy(baseline[, c(
  "id", "age_t1", "sex", "cannabis_current", "cannabis_frequency",
  "cudit_freq", "cudit_sev", "mmm_coping", "nsduh3"
), with = FALSE])

waves <- paste0("T", 12:17)
longitudinal <- data.table::rbindlist(lapply(seq_along(waves), function(wave_index) {
  wave <- waves[[wave_index]]
  wave_shift <- (wave_index - 1) * 0.10
  use_probability <- plogis(
    -0.30 + 0.65 * baseline$cannabis_current +
      0.25 * latent_externalizing - 0.05 * wave_shift
  )
  current <- rbinom(n, 1, use_probability)
  frequency <- ifelse(current == 1, sample(1:7, n, replace = TRUE), 0)
  cudit_c <- ifelse(
    current == 1,
    pmin(8, rpois(n, exp(0.65 + 0.08 * frequency + 0.15 * latent_externalizing))),
    0
  )
  coping <- ifelse(
    current == 1,
    pmin(20, pmax(0, rnorm(n, 6.5 + 0.9 * latent_adhd, 3))),
    0
  )
  cudit_p <- ifelse(
    current == 1,
    pmin(24, rpois(n, exp(0.25 + 0.15 * cudit_c + 0.035 * coping))),
    0
  )
  cud_probability <- plogis(
    -3.0 + 0.20 * cudit_c + 0.055 * coping +
      0.24 * latent_externalizing + 0.16 * latent_adhd +
      (0.05 + 0.06 * wave_index) * latent_somatic + random_intercept
  )
  cud <- ifelse(current == 1, rbinom(n, 1, cud_probability), 0)
  symptom_mean <- exp(
    -0.70 + 0.18 * cudit_c + 0.04 * coping +
      0.16 * latent_externalizing + 0.12 * latent_somatic + random_intercept / 2
  )
  symptoms <- ifelse(current == 1, pmin(11, rpois(n, symptom_mean)), 0)
  symptoms <- pmax(symptoms, ifelse(cud == 1, 2, 0))
  data.table::data.table(
    id = ids,
    wave = wave,
    age_wave = baseline$age_t1 + wave_index * 0.5,
    sex = baseline$sex,
    cannabis_current = current,
    cannabis_frequency = frequency,
    cudit_freq = cudit_c,
    cudit_sev = cudit_p,
    mmm_coping = coping,
    nsduh3 = sample(1:4, n, replace = TRUE),
    cud_case = cud,
    cud_symptom_count = symptoms,
    cud_symptom_n_observed = 11L
  )
}))

saveRDS(baseline, file.path(output_dir, "baseline_analysis.rds"))
saveRDS(
  baseline_longitudinal,
  file.path(output_dir, "baseline_longitudinal_predictors.rds")
)
saveRDS(longitudinal, file.path(output_dir, "longitudinal_analysis.rds"))
writeLines(
  c(
    paste0("Synthetic participants: ", n),
    paste0("Synthetic longitudinal rows: ", nrow(longitudinal)),
    "For software-contract testing only; never use for scientific interpretation."
  ),
  file.path(output_dir, "README.txt")
)
cat(normalizePath(output_dir), "\n")
