# PATH-CANN analysis v20: audited PGS selection for repeated CUD analyses.

select_cud_followup_pgs <- function(
    baseline_primary_results,
    pgs_manifest,
    output_dir,
    externalizing_pgs = "pgs_externalizing",
    cud_pgs = "pgs_cud") {
  check_required_columns(
    baseline_primary_results,
    c("pgs", "phenotype", "status", "passes_global_fdr"),
    "Baseline results used for CUD follow-up PGS selection"
  )
  check_required_columns(
    pgs_manifest, c("variable", "label", "domain"), "PGS manifest"
  )
  anchors <- c(externalizing_pgs, cud_pgs)
  missing_anchors <- setdiff(anchors, pgs_manifest$variable)
  if (length(missing_anchors)) {
    stop(
      "CUD follow-up anchor PGS(s) missing from the manifest: ",
      paste(missing_anchors, collapse = ", "), call. = FALSE
    )
  }

  hits <- unique(data.table::copy(baseline_primary_results[
    status == "ok" & passes_global_fdr %in% TRUE,
    .(pgs, phenotype)
  ]))
  externalizing_phenotypes <- unique(
    hits[pgs == externalizing_pgs, phenotype]
  )
  if (!length(externalizing_phenotypes)) {
    stop(
      "Externalizing must have at least one baseline global-FDR phenotype ",
      "to define the outside-externalizing selection rule.", call. = FALSE
    )
  }

  coverage <- hits[, .(
    baseline_fdr_phenotypes_n = data.table::uniqueN(phenotype),
    shared_with_externalizing_n = data.table::uniqueN(
      phenotype[phenotype %in% externalizing_phenotypes]
    ),
    outside_externalizing_n = data.table::uniqueN(
      phenotype[!phenotype %in% externalizing_phenotypes]
    )
  ), by = pgs]
  selection <- merge(
    data.table::copy(pgs_manifest)[, .(
      pgs = variable, pgs_label = label, pgs_domain = domain
    )],
    coverage, by = "pgs", all.x = TRUE
  )
  for (variable in c(
    "baseline_fdr_phenotypes_n",
    "shared_with_externalizing_n",
    "outside_externalizing_n"
  )) {
    data.table::set(selection, which(is.na(selection[[variable]])), variable, 0L)
  }
  selection[, selected_for_cud_followup :=
    pgs %in% anchors | outside_externalizing_n > 0]
  selection[, selection_reason := data.table::fcase(
    pgs == externalizing_pgs,
    "A priori externalizing anchor",
    pgs == cud_pgs,
    "A priori CUD comparator",
    outside_externalizing_n > 0,
    "At least one baseline FDR phenotype outside the externalizing hit set",
    default = paste0(
      "Not selected: no baseline FDR phenotype outside the externalizing hit set"
    )
  )]
  selection[, selection_order := data.table::fcase(
    pgs == externalizing_pgs, 1,
    pgs == cud_pgs, 2,
    default = 3
  )]
  data.table::setorder(
    selection, -selected_for_cud_followup, selection_order,
    -outside_externalizing_n, -baseline_fdr_phenotypes_n, pgs_label
  )
  selection[, selection_order := NULL]

  selected_ids <- selection[selected_for_cud_followup == TRUE, pgs]
  selected_manifest <- data.table::copy(
    pgs_manifest[match(selected_ids, variable)]
  )
  if (anyNA(selected_manifest$variable)) {
    stop("Selected CUD follow-up PGS manifest could not be resolved.", call. = FALSE)
  }
  write_csv(
    selection,
    file.path(output_dir, "longitudinal", "cud_followup_pgs_selection.csv")
  )
  list(audit = selection, manifest = selected_manifest)
}
