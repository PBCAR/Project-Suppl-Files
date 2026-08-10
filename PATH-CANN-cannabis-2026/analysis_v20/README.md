# PATH-CANN analysis v20

This is the corrected clean production package. Earlier analysis folders are unchanged.

## Run in RStudio

Open `run_all.R` and click **Source**. The script resolves the v20 folder from the active RStudio source file, creates a new additive output directory, and runs the complete production analysis.

The three input paths can be overridden before sourcing:

```r
Sys.setenv(PATH_CANN_BASELINE_FILE = "/full/path/baseline_analysis.rds")
Sys.setenv(PATH_CANN_BASELINE_LONGITUDINAL_FILE = "/full/path/baseline_longitudinal_predictors.rds")
Sys.setenv(PATH_CANN_LONGITUDINAL_FILE = "/full/path/longitudinal_analysis.rds")
```

The baseline input must contain the raw non-negative MPT parameters `Breakpoint`, `Intensity`, `Omax`, `Pmax`, and `Alpha`. v20 always recomputes the five analysis variables with `log1p`; it does not trust legacy logged columns.

## Locked analysis definitions

- Baseline screen: 12 PGSs by 30 cannabis phenotypes; one global BH family of 360 designated-population tests.
- Figure 1 displays all 30 phenotypes. The five MPT rows are explicitly labelled as current-user outcomes, and its coverage counts use the same complete phenotype set as Table 1 and PGS selection.
- The Figure 1 domain strip reads `Effects of cannabis use (current users)` because every phenotype in that domain is restricted to current users.
- Structural-zero all-participant outcomes: CUDIT-C, CUDIT-P, all five MMM motives, and expenditure. Their current-user estimates are effect-size sensitivities.
- Current-user outcomes: the five MPT parameters and cannabis-effect measures. Non-users are excluded from their primary analyses.
- Acceptability: medical and recreational `No opinion` responses are missing, leaving a 1–4 ordinal scale.
- Repeated CUD follow-up: externalizing and CUD anchors plus PGSs with a global-FDR phenotype outside the externalizing hit set.
- Wave-specific CUD estimates: marginal odds ratios at T12–T17 from a wave-saturated logistic model with participant-clustered HC1 covariance.
- Temporal heterogeneity: one omnibus PGS-by-wave Wald test per selected PGS; BH correction is separate from the 36 wave-specific tests.
- Figure S3: baseline genetic analysis cohort only at baseline and follow-up. A wave-level audit verifies its CUD denominator, cases, prevalence, and exclusion of non-genetic participants.

Every run writes `version.log`, `analysis_summary.txt`, coding/cohort audits, CSV results, and separate title-and-legend text files for figures and tables.
