# PATH-CANN analysis v20 figure module.
#
# This is an additive replacement candidate for R/05_figures.R. The original
# file remains unchanged. Plot titles and manuscript legends are written to
# separate text files so that the graphics themselves contain only essential
# axis, facet, and aesthetic labels.

save_analysis_plot <- function(plot, stem, width = 10, height = 7) {
  ggplot2::ggsave(
    paste0(stem, ".pdf"), plot,
    width = width, height = height, units = "in"
  )
  ggplot2::ggsave(
    paste0(stem, ".png"), plot,
    width = width, height = height, units = "in", dpi = 300
  )
}

write_display_text <- function(directory, stem, title, legend) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      paste0("Title: ", title),
      "",
      paste0("Legend: ", legend)
    ),
    file.path(directory, paste0(stem, "_title_and_legend.txt")),
    useBytes = TRUE
  )
}

baseline_figure_mpt_phenotypes <- c(
  "logBreakpoint", "logIntensity", "logOmax", "logPmax", "logAlpha"
)

baseline_fdr_hits <- function(primary_results) {
  unique(data.table::copy(primary_results[
    status == "ok" & passes_global_fdr %in% TRUE,
    .(
      pgs, pgs_label, pgs_domain,
      phenotype, phenotype_label, phenotype_domain
    )
  ]))
}

make_pgs_coverage_data <- function(primary_results) {
  hits <- baseline_fdr_hits(primary_results)
  if (!nrow(hits)) {
    stop("No baseline global-FDR hits were available for the coverage panel.", call. = FALSE)
  }
  externalizing_phenotypes <- unique(
    hits[pgs_label == "Externalizing", phenotype]
  )
  if (!length(externalizing_phenotypes)) {
    stop("Externalizing had no global-FDR phenotypes in the supplied results.", call. = FALSE)
  }

  coverage_counts <- hits[, .(
    total_fdr_phenotypes = data.table::uniqueN(phenotype),
    shared_with_externalizing = data.table::uniqueN(
      phenotype[phenotype %in% externalizing_phenotypes]
    ),
    outside_externalizing = data.table::uniqueN(
      phenotype[!phenotype %in% externalizing_phenotypes]
    )
  ), by = .(pgs, pgs_label, pgs_domain)]

  pgs_rows <- unique(data.table::copy(primary_results[, .(
    pgs, pgs_label, pgs_domain
  )]))
  coverage <- merge(
    pgs_rows, coverage_counts,
    by = c("pgs", "pgs_label", "pgs_domain"), all.x = TRUE
  )
  for (variable in c(
    "total_fdr_phenotypes",
    "shared_with_externalizing",
    "outside_externalizing"
  )) {
    data.table::set(coverage, which(is.na(coverage[[variable]])), variable, 0L)
  }

  coverage <- coverage[
    order(pgs_label != "Externalizing", -total_fdr_phenotypes, pgs_label)
  ]
  coverage[, pgs_order := seq_len(.N)]
  coverage
}

make_baseline_pgs_figure <- function(primary_results, output_dir) {
  figure_results <- data.table::copy(primary_results)
  plot_data <- data.table::copy(
    figure_results[status == "ok" & !is.na(p_value)]
  )
  if (!"passes_global_fdr" %in% names(plot_data)) {
    stop("Heatmap input must contain passes_global_fdr.", call. = FALSE)
  }

  plot_data[
    phenotype %in% baseline_figure_mpt_phenotypes,
    phenotype_label := paste0(phenotype_label, " (current users)")
  ]
  plot_data[
    phenotype_domain == "Effects of cannabis use",
    phenotype_domain := "Effects of cannabis use (current users)"
  ]
  phenotype_levels <- rev(unique(plot_data$phenotype_label))
  coverage <- make_pgs_coverage_data(figure_results)
  pgs_levels <- coverage$pgs_label
  plot_data[, phenotype_label := factor(phenotype_label, levels = phenotype_levels)]
  plot_data[, pgs_label := factor(pgs_label, levels = pgs_levels)]
  plot_data[, neg_log10_p := -log10(pmax(p_value, .Machine$double.xmin))]
  plot_data[, fdr_marker := data.table::fifelse(
    passes_global_fdr %in% TRUE, "*", ""
  )]

  heatmap_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = pgs_label, y = phenotype_label, fill = neg_log10_p)
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = fdr_marker),
      colour = "black", fontface = "bold", size = 4.2, vjust = 0.55
    ) +
    ggplot2::scale_fill_viridis_c(name = expression(-log[10](italic(p)))) +
    ggplot2::facet_grid(
      phenotype_domain ~ ., scales = "free_y", space = "free_y"
    ) +
    ggplot2::labs(x = "Polygenic score", y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank()
    )

  significant_pgs_levels <- coverage$pgs_label
  coverage_long <- data.table::rbindlist(list(
    coverage[, .(
      pgs_label,
      coverage_row = "Total FDR phenotypes",
      count = total_fdr_phenotypes
    )],
    coverage[, .(
      pgs_label,
      coverage_row = "Shared with externalizing",
      count = shared_with_externalizing
    )],
    coverage[, .(
      pgs_label,
      coverage_row = "Outside externalizing",
      count = outside_externalizing
    )]
  ))
  coverage_long[, pgs_label := factor(pgs_label, levels = significant_pgs_levels)]
  coverage_long[, coverage_row := factor(
    coverage_row,
    levels = c(
      "Outside externalizing",
      "Shared with externalizing",
      "Total FDR phenotypes"
    )
  )]
  coverage_long[, count_label := data.table::fifelse(
    count > 0, as.character(count), ""
  )]

  coverage_plot <- ggplot2::ggplot(
    coverage_long,
    ggplot2::aes(x = pgs_label, y = coverage_row)
  ) +
    ggplot2::geom_text(
      ggplot2::aes(label = count_label),
      size = 3.8, colour = "black", fontface = "bold"
    ) +
    ggplot2::labs(x = "Polygenic score", y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank()
    )

  combined <- patchwork::wrap_plots(
    heatmap_plot, coverage_plot,
    ncol = 1, heights = c(4.5, 1.15)
  ) +
    patchwork::plot_annotation(tag_levels = "A")

  stem <- "baseline_pgs_phenotype_heatmap"
  save_analysis_plot(
    combined,
    file.path(output_dir, "figures", stem),
    width = 14, height = 17
  )
  write_display_text(
    file.path(output_dir, "figures"),
    stem,
    "Figure 1. Polygenic association coverage across 30 baseline cannabis-related phenotypes.",
    paste0(
      "Panel A shows -log10 p-values for ", nrow(plot_data), " associations involving ",
      data.table::uniqueN(plot_data$pgs), " PGSs and ",
      data.table::uniqueN(plot_data$phenotype), " phenotypes. The five Marijuana Purchase Task ",
      "phenotypes are labelled as current-user outcomes and use log1p-transformed raw parameters; ",
      "non-users are not imputed for these measures. The Effects of cannabis use domain is also ",
      "labelled as current users because every phenotype in that domain is user-defined. An ",
      "asterisk marks an association passing the ",
      "single global Benjamini-Hochberg FDR at q < 0.05 across the full prespecified ",
      unique(primary_results$fdr_family_tests), "-test family. ",
      "Panel B summarizes for all 12 PGSs the total number of ",
      "associated phenotypes, the number also associated with externalizing, and the number outside ",
      "externalizing's complete 30-phenotype hit set. Positive counts are printed directly; zero ",
      "cells are left blank. ",
      "Counts describe association coverage and do not establish ",
      "independent or unique genetic effects."
    )
  )
  combined
}

make_liability_evidence_table <- function(
    primary_results, cud_pgs_results, output_dir) {
  pgs_rows <- unique(data.table::copy(primary_results[, .(
    pgs, pgs_label, pgs_domain
  )]))
  hits <- baseline_fdr_hits(primary_results)
  coverage <- make_pgs_coverage_data(primary_results)[, .(
    pgs,
    baseline_fdr_phenotypes_n = total_fdr_phenotypes,
    shared_with_externalizing_n = shared_with_externalizing,
    outside_externalizing_n = outside_externalizing
  )]
  evidence <- merge(pgs_rows, coverage, by = "pgs", all.x = TRUE)
  for (variable in c(
    "baseline_fdr_phenotypes_n",
    "shared_with_externalizing_n",
    "outside_externalizing_n"
  )) {
    data.table::set(evidence, which(is.na(evidence[[variable]])), variable, 0L)
  }
  evidence[, baseline_any_fdr := baseline_fdr_phenotypes_n > 0]
  evidence[, current_cannabis_use_fdr := pgs %in%
    hits[phenotype == "cannabis_current", pgs]]
  evidence[, cudit_p_fdr := pgs %in%
    hits[phenotype == "cudit_sev", pgs]]

  cud_results <- data.table::copy(cud_pgs_results$wave_effects)
  interaction_results <- data.table::copy(cud_pgs_results$interactions)
  selection <- data.table::copy(cud_pgs_results$selection)
  required_cud_columns <- c(
    "pgs", "wave", "passes_wave_fdr", "fdr_family_tests"
  )
  missing_cud_columns <- setdiff(required_cud_columns, names(cud_results))
  if (length(missing_cud_columns)) {
    stop(
      "CUD-status results are missing: ",
      paste(missing_cud_columns, collapse = ", "),
      call. = FALSE
    )
  }
  check_required_columns(
    selection,
    c("pgs", "selected_for_cud_followup", "selection_reason"),
    "CUD follow-up PGS selection audit"
  )
  evidence <- merge(
    evidence,
    selection[, .(
      pgs, selected_for_cud_followup, cud_followup_selection_reason = selection_reason
    )],
    by = "pgs", all.x = TRUE
  )
  cud_hits <- cud_results[passes_wave_fdr %in% TRUE]
  later_summary <- cud_hits[, .(
    later_cud_fdr_waves_n = data.table::uniqueN(wave),
    later_cud_fdr_waves = paste(
      unique(wave)[order(as.integer(sub("^T", "", unique(wave))))],
      collapse = ", "
    )
  ), by = pgs]
  evidence <- merge(evidence, later_summary, by = "pgs", all.x = TRUE)
  evidence <- merge(
    evidence,
    interaction_results[, .(
      pgs,
      pgs_by_wave_p_value = p_value,
      pgs_by_wave_fdr_q = fdr_q_interaction,
      pgs_by_wave_passes_fdr = passes_interaction_fdr
    )],
    by = "pgs", all.x = TRUE
  )
  evidence[
    selected_for_cud_followup == TRUE & is.na(later_cud_fdr_waves_n),
    later_cud_fdr_waves_n := 0L
  ]
  evidence[
    selected_for_cud_followup == TRUE & is.na(later_cud_fdr_waves),
    later_cud_fdr_waves := "None"
  ]
  evidence[
    selected_for_cud_followup == FALSE,
    `:=`(later_cud_fdr_waves_n = NA_integer_, later_cud_fdr_waves = "Not evaluated")
  ]
  evidence[, later_cud_any_fdr := data.table::fifelse(
    selected_for_cud_followup,
    later_cud_fdr_waves_n > 0,
    NA
  )]
  evidence <- evidence[
    order(
      pgs != "pgs_externalizing",
      pgs != "pgs_cud",
      -selected_for_cud_followup,
      -baseline_fdr_phenotypes_n,
      -later_cud_fdr_waves_n,
      pgs_label
    )
  ]
  evidence[, polygenic_score := pgs_label]
  table_output <- evidence[, .(
    polygenic_score,
    pgs_domain,
    baseline_fdr_phenotypes_n,
    shared_with_externalizing_n,
    outside_externalizing_n,
    baseline_any_fdr,
    current_cannabis_use_fdr,
    cudit_p_fdr,
    selected_for_cud_followup,
    cud_followup_selection_reason,
    pgs_by_wave_p_value,
    pgs_by_wave_fdr_q,
    pgs_by_wave_passes_fdr,
    later_cud_fdr_waves_n,
    later_cud_fdr_waves,
    later_cud_any_fdr
  )]

  table_dir <- file.path(output_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    table_output,
    file.path(table_dir, "table_1_polygenic_liability_evidence.csv"),
    na = ""
  )
  write_display_text(
    table_dir,
    "table_1_polygenic_liability_evidence",
    "Table 1. Polygenic association coverage across baseline cannabis phenotypes and later CUD status.",
    paste0(
      "Baseline counts are distinct PGS-phenotype associations passing the global BH FDR in the ",
      "360-test family. Shared-with-externalizing counts phenotype types also associated with the ",
      "externalizing PGS; outside-externalizing counts phenotype types absent from the externalizing ",
      "hit set. Externalizing and the CUD PGS were retained a priori for follow-up; another PGS was ",
      "retained when it had at least one baseline FDR phenotype outside the externalizing hit set. ",
      "Later CUD results come from marginal logistic models across T12-T17 with participant-clustered ",
      "HC1 covariance. Omnibus PGS-by-wave Wald tests are BH-corrected across the ",
      unique(interaction_results$fdr_family_tests), " selected scores; wave-specific estimates are ",
      "BH-corrected across ", unique(cud_results$fdr_family_tests), " tests. In total, ",
      nrow(table_output), " PGSs were tested at baseline, ",
      sum(table_output$baseline_any_fdr), " were associated with at least one baseline phenotype, and ",
      sum(table_output$selected_for_cud_followup), " were retained for CUD follow-up. Among retained ",
      "scores, ", sum(table_output$later_cud_any_fdr, na.rm = TRUE),
      " were associated with CUD status at one or more waves. ",
      "Coverage counts do not establish mutually independent polygenic effects."
    )
  )
  table_output
}

make_global_hit_population_comparison_figure <- function(comparison, output_dir) {
  if (!nrow(comparison)) return(NULL)
  plot_data <- data.table::copy(comparison[status == "ok"])
  if (!nrow(plot_data)) return(NULL)

  plot_data[, population_label := factor(
    comparison_population,
    levels = c("All participants", "Current users only")
  )]
  plot_data[, association := paste(pgs_label, "to", phenotype_label)]
  panel_labels <- c(
    ordinal = "A) Ordinal outcomes: proportional odds ratio",
    count = "B) Count outcomes: rate ratio",
    continuous = "C) Continuous outcomes: standardized coefficient"
  )
  unexpected_families <- setdiff(unique(plot_data$model_family), names(panel_labels))
  if (length(unexpected_families)) {
    stop(
      "Unsupported comparison model family/families: ",
      paste(unexpected_families, collapse = ", "),
      call. = FALSE
    )
  }
  plot_data[, effect_panel := unname(panel_labels[model_family])]
  plot_data[, `:=`(
    plot_estimate = data.table::fifelse(
      model_family == "continuous", as.numeric(estimate), as.numeric(effect_ratio)
    ),
    plot_low = data.table::fifelse(
      model_family == "continuous", as.numeric(conf_low), as.numeric(ratio_low)
    ),
    plot_high = data.table::fifelse(
      model_family == "continuous", as.numeric(conf_high), as.numeric(ratio_high)
    )
  )]
  plot_data <- plot_data[
    is.finite(plot_estimate) & is.finite(plot_low) & is.finite(plot_high)
  ]
  plot_data[, effect_panel := factor(effect_panel, levels = unname(panel_labels))]
  association_order <- unique(
    plot_data[order(effect_panel, selection_p_value)]$association
  )
  plot_data[, association := factor(association, levels = rev(association_order))]
  references <- unique(plot_data[, .(effect_panel)])
  references[, reference := data.table::fifelse(
    effect_panel == panel_labels[["continuous"]], 0, 1
  )]
  dodge <- ggplot2::position_dodge(width = 0.55)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = plot_estimate, y = association,
      colour = population_label, shape = population_label
    )
  ) +
    ggplot2::geom_vline(
      data = references,
      ggplot2::aes(xintercept = reference),
      linetype = 3, colour = "grey55", inherit.aes = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = plot_low, xmax = plot_high),
      position = dodge, orientation = "y", width = 0.18, linewidth = 0.6
    ) +
    ggplot2::geom_point(position = dodge, size = 2.5) +
    ggplot2::facet_wrap(
      ggplot2::vars(effect_panel), nrow = 1, scales = "free"
    ) +
    ggplot2::labs(
      x = "Effect per 1-SD higher PGS", y = NULL,
      colour = NULL, shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10.5) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.spacing.x = grid::unit(1.2, "lines"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank()
    )

  stem <- "global_fdr_hits_all_vs_users_effects"
  save_analysis_plot(
    plot, file.path(output_dir, "figures", stem), width = 20, height = 10
  )
  write_display_text(
    file.path(output_dir, "figures"),
    stem,
    "Figure S1. Global-FDR PGS-phenotype associations in all participants and current users.",
    paste0(
      "Panel A shows proportional odds ratios for ordinal outcomes, Panel B rate ratios for count ",
      "outcomes, and Panel C standardized coefficients for continuous outcomes. Points and bars are ",
      "estimates and 95% confidence intervals. Panels are separated because their effect measures are ",
      "not directly comparable. Cannabis-use-status outcomes are excluded because they are constant ",
      "among current users. Marijuana Purchase Task outcomes are excluded because they are defined ",
      "only among current users, so the two population labels would not represent different analysis ",
      "samples. The non-primary population estimate is an effect-size-only sensitivity."
    )
  )
  plot
}

make_externalizing_retention_figure <- function(comparison, output_dir) {
  zero_coded_phenotypes <- c(
    "cudit_freq", "cudit_sev", "cannabis_expenditure",
    "mmm_enhancement", "mmm_conformity", "mmm_expansion",
    "mmm_coping", "mmm_social"
  )
  externalizing <- data.table::copy(comparison[
    status == "ok" &
      pgs_label == "Externalizing" &
      phenotype %in% zero_coded_phenotypes
  ])
  if (!nrow(externalizing)) {
    stop("No externalizing global-FDR zero-coded comparison pairs were available.", call. = FALSE)
  }

  all_rows <- externalizing[
    comparison_population == "All participants",
    .(
      phenotype, phenotype_label, model_family,
      estimate_all = as.numeric(estimate),
      conf_low_all = as.numeric(conf_low),
      conf_high_all = as.numeric(conf_high)
    )
  ]
  user_rows <- externalizing[
    comparison_population == "Current users only",
    .(
      phenotype, phenotype_label, model_family,
      estimate_users = as.numeric(estimate),
      conf_low_users = as.numeric(conf_low),
      conf_high_users = as.numeric(conf_high)
    )
  ]
  paired <- merge(
    all_rows, user_rows,
    by = c("phenotype", "phenotype_label", "model_family"),
    all = FALSE
  )
  paired <- paired[
    is.finite(estimate_all) & abs(estimate_all) > 1e-12 &
      is.finite(estimate_users) &
      is.finite(conf_low_all) & is.finite(conf_high_all) &
      is.finite(conf_low_users) & is.finite(conf_high_users)
  ]
  if (!nrow(paired)) {
    stop("Externalizing retention estimates could not be normalized.", call. = FALSE)
  }

  paired[, `:=`(
    estimate_retained_percent = 100 * estimate_users / estimate_all,
    all_ci_a = 100 * conf_low_all / estimate_all,
    all_ci_b = 100 * conf_high_all / estimate_all,
    users_ci_a = 100 * conf_low_users / estimate_all,
    users_ci_b = 100 * conf_high_users / estimate_all
  )]
  paired[, `:=`(
    all_plot_low = pmin(all_ci_a, all_ci_b),
    all_plot_high = pmax(all_ci_a, all_ci_b),
    users_plot_low = pmin(users_ci_a, users_ci_b),
    users_plot_high = pmax(users_ci_a, users_ci_b),
    user_ci_includes_null = conf_low_users <= 0 & conf_high_users >= 0,
    direction_consistent = sign(estimate_all) == sign(estimate_users)
  )]
  paired[, user_ci_overlap_length := pmax(
    0,
    pmin(conf_high_all, conf_high_users) - pmax(conf_low_all, conf_low_users)
  )]
  paired[, user_ci_width := conf_high_users - conf_low_users]
  paired[, fraction_user_ci_overlapping_all_ci := data.table::fifelse(
    user_ci_width > 0,
    pmin(1, user_ci_overlap_length / user_ci_width),
    NA_real_
  )]

  phenotype_order <- paired[order(estimate_retained_percent)]$phenotype_label
  paired[, phenotype_label := factor(
    phenotype_label, levels = rev(phenotype_order)
  )]
  plot_data <- data.table::rbindlist(list(
    paired[, .(
      phenotype_label,
      population = "All participants",
      plot_estimate = 100,
      plot_low = all_plot_low,
      plot_high = all_plot_high
    )],
    paired[, .(
      phenotype_label,
      population = "Current users only",
      plot_estimate = estimate_retained_percent,
      plot_low = users_plot_low,
      plot_high = users_plot_high
    )]
  ))
  plot_data[, population := factor(
    population, levels = c("All participants", "Current users only")
  )]
  dodge <- ggplot2::position_dodge(width = 0.55)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = plot_estimate, y = phenotype_label,
      colour = population, shape = population
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.45) +
    ggplot2::geom_vline(xintercept = 50, colour = "grey55", linetype = 3) +
    ggplot2::geom_vline(xintercept = 100, colour = "grey35", linetype = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = plot_low, xmax = plot_high),
      orientation = "y", width = 0.16, linewidth = 0.65,
      position = dodge
    ) +
    ggplot2::geom_point(size = 2.7, position = dodge) +
    ggplot2::scale_colour_manual(values = c(
      "All participants" = "#E76F51",
      "Current users only" = "#2A9D8F"
    )) +
    ggplot2::labs(
      x = "Coefficient relative to all-participant estimate (%)",
      y = NULL, colour = NULL, shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank()
    )

  stem <- "externalizing_zero_coded_effect_retention"
  save_analysis_plot(
    plot, file.path(output_dir, "figures", stem), width = 10, height = 6.5
  )
  write_display_text(
    file.path(output_dir, "figures"),
    stem,
    "Figure 2. Retention of externalizing associations among current cannabis users.",
    paste0(
      "The figure includes baseline externalizing associations passing global FDR for phenotypes ",
      "structurally zero-coded among confirmed non-users. Estimates and 95% confidence intervals are ",
      "expressed relative to the corresponding all-participant coefficient, which is set to 100%. ",
      "The dashed lines mark 50% and 100% retention. Normalized intervals describe uncertainty in each ",
      "fitted coefficient while treating the all-participant point estimate as the scaling reference; ",
      "they are descriptive and do not test the difference between the overlapping samples."
    )
  )

  audit <- paired[, .(
    phenotype,
    phenotype_label = as.character(phenotype_label),
    model_family,
    estimate_all,
    conf_low_all,
    conf_high_all,
    estimate_users,
    conf_low_users,
    conf_high_users,
    estimate_retained_percent,
    user_ci_includes_null,
    direction_consistent,
    fraction_user_ci_overlapping_all_ci
  )]
  table_dir <- file.path(output_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    audit,
    file.path(table_dir, "externalizing_zero_coded_retention.csv"),
    na = ""
  )
  write_display_text(
    table_dir,
    "externalizing_zero_coded_retention",
    "Supplementary table. Externalizing coefficient retention and descriptive CI overlap among current users.",
    paste0(
      "Retention is 100 multiplied by the current-user coefficient divided by the all-participant ",
      "coefficient on the fitted model scale. The CI-overlap fraction is the length of intersection ",
      "between the two coefficient intervals divided by the width of the current-user interval. It is ",
      "reported descriptively only: because current users are nested within the all-participant sample, ",
      "CI overlap is not a test of effect equality and is strongly influenced by user-only precision."
    )
  )
  list(plot = plot, data = audit)
}

validate_longitudinal_profile <- function(profile) {
  prevalence <- data.table::copy(profile[
    wave != "T1" &
      measure %in% c("cannabis_current", "cud_case_plot") &
      population == "all"
  ])
  if (!nrow(prevalence)) {
    stop("No follow-up prevalence rows were available.", call. = FALSE)
  }
  prevalence[, signature := paste(
    n, signif(estimate, 15), signif(conf_low, 15), signif(conf_high, 15),
    sep = "|"
  )]
  check <- prevalence[, .(
    wave_count = data.table::uniqueN(wave),
    summary_count = data.table::uniqueN(signature)
  ), by = measure]
  failed <- check[wave_count > 1 & summary_count == 1]
  if (nrow(failed)) {
    stop(
      "Identical follow-up summaries were detected for ",
      paste(failed$measure, collapse = ", "),
      ". Rebuild the longitudinal profile after correcting wave-specific subsetting.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

make_longitudinal_profile_figure <- function(profile, output_dir) {
  validate_longitudinal_profile(profile)
  wave_levels <- unique(profile[order(wave_number)]$wave)
  bar_data <- data.table::copy(profile[
    measure %in% c("cannabis_current", "cud_case_plot") & population == "all"
  ])
  bar_data[, wave := factor(wave, levels = wave_levels)]
  bar_data[, prevalence_measure := factor(
    measure_label, levels = c("Current cannabis use", "DSM-5 CUD")
  )]
  bar_data[, label_y := pmin(0.98, conf_high + 0.035)]
  bar_ceiling <- min(
    1,
    max(0.50, max(bar_data$label_y, na.rm = TRUE) + 0.04)
  )
  prevalence_dodge <- ggplot2::position_dodge(width = 0.72)
  bar_plot <- ggplot2::ggplot(
    bar_data,
    ggplot2::aes(x = wave, y = estimate, fill = prevalence_measure)
  ) +
    ggplot2::geom_col(
      width = 0.68, position = prevalence_dodge, na.rm = TRUE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = conf_low, ymax = conf_high),
      width = 0.15, linewidth = 0.6,
      position = prevalence_dodge, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = label_y, label = sprintf("%.1f%%", 100 * estimate)),
      vjust = 0, size = 3.1, position = prevalence_dodge, na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = c("Current cannabis use" = "#2A9D8F", "DSM-5 CUD" = "#E76F51"),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, bar_ceiling),
      labels = function(x) paste0(round(100 * x), "%"),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(x = NULL, y = "Prevalence", fill = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank()
    )

  line_data <- data.table::copy(profile[
    (measure %in% c("cudit_freq", "cudit_sev") & population == "current_users") |
      (measure == "cud_symptom_plot" &
         population == "current_users_complete_symptoms")
  ])
  line_data[, wave := factor(wave, levels = wave_levels)]
  line_data[, measure_label := factor(
    measure_label,
    levels = c("CUDIT-C", "CUDIT-P", "DSM-5 CUD symptoms")
  )]
  line_data[, `:=`(
    plot_estimate = estimate,
    plot_low = conf_low,
    plot_high = conf_high
  )]
  line_ceiling <- max(5, 1.10 * max(line_data$plot_high, na.rm = TRUE))
  line_plot <- ggplot2::ggplot(
    line_data,
    ggplot2::aes(
      x = wave, y = plot_estimate,
      colour = measure_label, shape = measure_label,
      linetype = measure_label, group = measure_label
    )
  ) +
    ggplot2::geom_line(linewidth = 0.85, na.rm = TRUE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = plot_low, ymax = plot_high),
      width = 0.10, linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2.8, stroke = 0.9, na.rm = TRUE) +
    ggplot2::scale_colour_manual(
      values = c(
        "CUDIT-C" = "#0072B2",
        "CUDIT-P" = "#D55E00",
        "DSM-5 CUD symptoms" = "#7B3294"
      ), drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = c("CUDIT-C" = 16, "CUDIT-P" = 17, "DSM-5 CUD symptoms" = 15),
      drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "CUDIT-C" = "solid",
        "CUDIT-P" = "dashed",
        "DSM-5 CUD symptoms" = "dotdash"
      ), drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, line_ceiling),
      breaks = scales::pretty_breaks(n = 6),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(
      x = "Assessment wave", y = "Mean score",
      colour = NULL, shape = NULL, linetype = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank()
    )

  combined <- patchwork::wrap_plots(
    bar_plot, line_plot, ncol = 1, heights = c(1, 1.25)
  ) +
    patchwork::plot_annotation(tag_levels = "A")
  stem <- "longitudinal_cannabis_profile_current_users"
  save_analysis_plot(
    combined, file.path(output_dir, "figures", stem), width = 11, height = 8.5
  )
  write_display_text(
    file.path(output_dir, "figures"),
    stem,
    "Figure S3. Cannabis use, CUD prevalence and CUD-related measures across assessments in the genetic analysis cohort.",
    paste0(
      "The sample is restricted to participants in the baseline genetic analysis cohort. Panel A ",
      "shows current cannabis use and DSM-5 CUD prevalence with 95% confidence intervals; ",
      "labels are positioned above the upper interval. Panel B shows mean CUDIT-C, CUDIT-P and DSM-5 ",
      "CUD symptom counts among current users using distinct colours, shapes and line types. The three ",
      "raw scores have different theoretical ranges (0-8, 0-24 and 0-11) and their vertical levels are ",
      "therefore not directly comparable. DSM-5 symptom means require a complete 11-item assessment."
    )
  )
  combined
}

make_marginal_wave_cud_figure <- function(
    cud_results, interaction_results, output_dir) {
  plot_data <- data.table::copy(cud_results[
    status == "ok" &
      is.finite(effect_ratio) & is.finite(ratio_low) &
      is.finite(ratio_high) & ratio_low > 0 & ratio_high > 0
  ])
  if (!nrow(plot_data)) {
    stop("No finite marginal wave-specific CUD estimates were available.", call. = FALSE)
  }
  wave_levels <- unique(plot_data[order(as.integer(sub("^T", "", wave)))]$wave)
  pgs_levels <- unique(plot_data$pgs_label)
  plot_data[, `:=`(
    pgs_label = factor(pgs_label, levels = pgs_levels),
    wave = factor(wave, levels = wave_levels),
    fdr_marker = data.table::fifelse(
      passes_wave_fdr %in% TRUE, "*", ""
    )
  )]

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = wave, y = effect_ratio)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = 3, colour = "grey55") +
    ggplot2::geom_line(
      ggplot2::aes(group = 1), linewidth = 0.65, colour = "#145A64"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ratio_low, ymax = ratio_high),
      width = 0.12, linewidth = 0.65, colour = "#145A64"
    ) +
    ggplot2::geom_point(size = 2.8, colour = "#145A64") +
    ggplot2::geom_text(
      ggplot2::aes(y = ratio_high * 1.10, label = fdr_marker),
      colour = "black", size = 4.5, fontface = "bold"
    ) +
    ggplot2::scale_y_log10(
      expand = ggplot2::expansion(mult = c(0.08, 0.20))
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(pgs_label), ncol = 3
    ) +
    ggplot2::labs(
      x = "CUD assessment wave",
      y = "Marginal odds ratio per 1-SD higher PGS (log scale)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank()
    )

  stem <- "selected_pgs_wave_specific_cud_status"
  save_analysis_plot(
    plot, file.path(output_dir, "figures", stem), width = 12, height = 8.5
  )
  write_display_text(
    file.path(output_dir, "figures"),
    stem,
    "Figure 3. Wave-specific associations of selected polygenic scores with DSM-5 CUD status.",
    paste0(
      "Points and bars are adjusted population-level odds ratios and participant-clustered 95% ",
      "confidence intervals per 1-SD higher PGS at T12-T17. Each selected PGS is estimated in a ",
      "pooled logistic model saturated by wave for the PGS and every adjustment covariate; robust ",
      "HC1 standard errors account for repeated observations within participants. Within each PGS ",
      "panel, lines connect the repeated cross-sectional estimates to guide comparison across waves; ",
      "they do not represent individual trajectories. Externalizing and the CUD PGS were retained a ",
      "priori; another PGS ",
      "was retained when it had at least one baseline global-FDR phenotype outside the ",
      "externalizing-associated phenotype set. Asterisks identify wave-specific estimates passing ",
      "BH FDR across ", unique(plot_data$fdr_family_tests), " selected-PGS-by-wave tests. Omnibus ",
      "PGS-by-wave Wald tests use the joint participant-clustered covariance and are BH-corrected ",
      "across ", unique(interaction_results$fdr_family_tests), " selected PGSs; they are reported ",
      "separately. Differences in pointwise significance alone do not demonstrate heterogeneity."
    )
  )
  plot
}

run_figures <- function(
    baseline_results,
    cud_pgs_results,
    longitudinal_profile,
    output_dir) {
  evidence_table <- make_liability_evidence_table(
    baseline_results$primary, cud_pgs_results, output_dir
  )
  retention <- make_externalizing_retention_figure(
    baseline_results$population_comparison, output_dir
  )
  list(
    baseline = make_baseline_pgs_figure(
      baseline_results$primary, output_dir
    ),
    externalizing_retention = retention$plot,
    liability_evidence_table = evidence_table,
    population_comparison = make_global_hit_population_comparison_figure(
      baseline_results$population_comparison, output_dir
    ),
    longitudinal_profile = make_longitudinal_profile_figure(
      longitudinal_profile, output_dir
    ),
    wave_specific_cud = make_marginal_wave_cud_figure(
      cud_pgs_results$wave_effects,
      cud_pgs_results$interactions,
      output_dir
    )
  )
}
