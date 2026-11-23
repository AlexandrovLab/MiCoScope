# Reproducible execution and diagnostic reporting, not a new correction algorithm.
configure_batch_rng <- function(seed, threads) {
  if (length(seed) != 1L || is.na(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed))
    stop("seed must be an integer between 0 and 2147483647")
  if (length(threads) != 1L || is.na(threads) || threads < 1 || threads != floor(threads))
    stop("threads must be a positive integer")
  if (.Platform$OS.type != "unix" && threads > 1L)
    stop("Reproducible parallel ConQuR currently requires the Unix multicore backend")
  RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(as.integer(seed))
  invisible(TRUE)
}

quality_summary <- function(baseline, final) {
  ok <- length(baseline) == 1L && length(final) == 1L && is.finite(baseline) && is.finite(final)
  ratio <- if (ok && baseline != 0) (baseline - final) / baseline else NA_real_
  label <- if (!ok) "not_evaluable_nonfinite_R2" else if (baseline <= 0)
    "not_interpretable_nonpositive_baseline" else if (final < 0)
    "not_interpretable_negative_final_R2" else if (final < baseline)
    "lower_batch_R2_descriptive_only" else "no_reduction_in_batch_R2"
  interpretable <- ok && baseline > 0 && final >= 0
  data.frame(baseline_R2 = if (length(baseline) == 1) baseline else NA_real_,
             final_R2 = if (length(final) == 1) final else NA_real_,
             diagnostic_relative_change = ratio,
             reported_reduction_percent = if (interpretable) ratio * 100 else NA_real_,
             quality_status = label, biological_signal_preservation = "not_assessed",
             taxon_fit_diagnostics = "not_exposed_by_pinned_dependency", stringsAsFactors = FALSE)
}

report_reduction_percent <- function(baseline, final) quality_summary(baseline, final)$reported_reduction_percent

# Which table to keep after correction. Never stops the run and never selects on
# NA/NaN/Inf: a candidate whose batch R2 could not be evaluated is not selectable; a non-positive or
# unevaluable baseline means there is no measurable batch effect to remove; when no evaluable
# candidate improves on the baseline, the uncorrected table is kept. Every outcome is recorded.
#   baseline: batch R2 before correction; candidates: named numeric batch R2 per candidate, in
#   preference order (an earlier candidate wins ties). threshold applies to the first candidate only
#   (the existing stop-at-phase-1 rule); it is reported, not used as an acceptance criterion.
#   baseline_p: p of the baseline batch term on relative-abundance Bray; when it is
#   finite and >= alpha there is no measurable batch effect and the uncorrected table is kept.
decide_correction <- function(baseline, candidates, threshold = 0.25, baseline_p = NA_real_, alpha = 0.05) {
  finite <- function(x) length(x) == 1L && !is.na(x) && is.finite(x)
  out <- function(choice, decision, reason) list(choice = choice, decision = decision, reason = reason)
  if (!finite(baseline))
    return(out(NA_character_, "uncorrected_not_evaluable", "baseline batch R2 could not be evaluated"))
  if (baseline <= 0)
    return(out(NA_character_, "uncorrected_no_batch_effect",
               sprintf("no measurable batch effect before correction (baseline batch R2 = %g)", baseline)))
  if (finite(baseline_p) && baseline_p >= alpha)
    return(out(NA_character_, "uncorrected_no_batch_effect",
               sprintf("batch effect not significant before correction (relative-abundance Bray R2 = %g, p = %g >= %g)",
                       baseline, baseline_p, alpha)))
  ok <- vapply(candidates, finite, logical(1))
  skipped <- names(candidates)[!ok]
  note <- if (length(skipped)) paste0("; not evaluable: ", paste(skipped, collapse = ", ")) else ""
  usable <- candidates[ok & candidates < baseline]
  if (!length(usable))
    return(out(NA_character_, "uncorrected_no_improvement",
               paste0("no evaluable candidate lowered batch R2 below the baseline ", sprintf("%g", baseline), note)))
  best <- names(usable)[which.min(usable)]
  out(best, "corrected", sprintf("%s lowered batch R2 from %g to %g (absolute change %g)%s",
                                 best, baseline, usable[[best]], usable[[best]] - baseline, note))
}

# Effect of correction on the sample-type comparison. Per taxon, log2 ratio of the
# mean relative abundance in group A over group B, in the input table and in the table used downstream.
# Computed only when exactly two sample types remain; otherwise one row says why. Group A is the declared
# tumour value when it is one of the two, else the first in sorted order.
# protected_covariate is NA on the sample-type rows and yes/no on the rows of a user-named comparison, which
# are appended to the same file (write_named_effect_sizes).
write_effect_sizes <- function(tables, meta, group_col, tumour_values, path, append = FALSE, protected = NA,
                               reason = NULL, what = "sample types") {
  after_tag <- intersect(c("ConQuR_tuned", "UNCORRECTED_by_decision", "UNCORRECTED_passthrough"), names(tables))[1]
  groups <- if (!is.null(group_col) && nzchar(group_col) && group_col %in% names(meta))
    sort(unique(trimws(as.character(meta[[group_col]])))) else character(0)
  put <- function(df) write.table(cbind(df, protected_covariate = protected), path, sep = "\t", quote = FALSE,
                                  row.names = FALSE, na = "NA", append = append, col.names = !append)
  empty <- function(reason) put(data.frame(taxon = NA, group_column = if (is.null(group_col)) NA else group_col,
    group_a = NA, group_b = NA, log2_ratio_before = NA, log2_ratio_after = NA, status = paste("not computed:", reason)))
  if (!is.null(reason)) return(invisible(empty(reason)))
  if (is.na(after_tag)) return(invisible(empty("no table used downstream")))
  if (length(groups) != 2L) return(invisible(empty(sprintf("needs exactly two %s (found %d)", what, length(groups)))))
  a <- if (any(toupper(groups) %in% toupper(tumour_values))) groups[toupper(groups) %in% toupper(tumour_values)][1] else groups[1]
  b <- setdiff(groups, a)
  ratio <- function(x) {
    x <- as.matrix(x); g <- trimws(as.character(meta[rownames(x), group_col])); rel <- x / rowSums(x)
    ma <- colMeans(rel[g %in% a, , drop = FALSE]); mb <- colMeans(rel[g %in% b, , drop = FALSE])
    ifelse(ma > 0 & mb > 0, log2(ma / mb), NA_real_)
  }
  before <- ratio(tables$raw); after <- ratio(tables[[after_tag]])[names(before)]
  put(data.frame(taxon = names(before), group_column = group_col, group_a = a, group_b = b,
                 log2_ratio_before = before, log2_ratio_after = after, status = "computed"))
  invisible(TRUE)
}

# The same before/after effect sizes for a user-named metadata column (--effect_column), appended to the
# sample-type rows. The column must exist and hold exactly two non-missing values among the analysis
# samples; otherwise one row says why. Samples with a missing value are left out of this comparison only.
# The column is protected when it is one of the covariates used for correction.
write_named_effect_sizes <- function(tables, meta, col, covariates, type_column, path) {
  protected <- if (col %in% covariates) "yes" else "no"
  go <- function(...) write_effect_sizes(tables, meta, col, character(0), path, append = TRUE, protected = protected, ...)
  if (!col %in% names(meta)) return(go(reason = "column not found in the metadata"))
  if (identical(col, type_column)) return(go(reason = "same column as the sample-type comparison"))
  v <- trimws(as.character(meta[[col]]))
  meta <- meta[!(is.na(v) | v %in% c("", "NA")), , drop = FALSE]
  go(what = "non-missing values among the analysis samples")
}

write_batch_records <- function(result, tab, predictors, args) {
  out <- dirname(args$prefix)
  q <- quality_summary(result$baseline_r2, result$final_r2)
  q$execution_status <- if (!is.null(result$decision) && grepl("^uncorrected_", result$decision)) result$decision else
    if (isTRUE(result$success)) "completed" else "uncorrected_explicit_fallback"
  q$absolute_change <- if (length(result$baseline_r2) == 1L && length(result$final_r2) == 1L)
    as.numeric(result$final_r2) - as.numeric(result$baseline_r2) else NA_real_
  q$correction_decision_reason <- if (is.null(result$decision_reason)) NA_character_ else result$decision_reason
  # p of the baseline batch term on relative-abundance Bray, the value the correction decision used.
  db <- if (exists("decision_baseline", envir = .GlobalEnv)) get("decision_baseline", envir = .GlobalEnv) else NULL
  q$baseline_p <- if (is.null(db) || is.null(db$batch_p)) NA_real_ else as.numeric(db$batch_p)
  q$sample_type_note <- if (is.null(result$sample_type_note)) NA_character_ else result$sample_type_note
  q$confounding_note <- if (is.null(result$confounding_note) || is.na(result$confounding_note)) NA_character_ else result$confounding_note
  q$n_samples <- nrow(tab)
  q$n_taxa <- ncol(tab)
  q$method <- args$method
  write.table(q, file.path(out, "batch_correction_quality.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  controls <- if (is.null(args$control_values) || !nzchar(args$control_values) || toupper(args$control_values) == "NULL") "none declared" else
    paste0(if (isTRUE(args$include_controls)) "included: " else "excluded: ", args$control_values)
  settings <- data.frame(key = c("seed", "threads", "rng_kind", "normal_kind", "sample_kind", "covariates", "batch_column", "method", "comparison_scope", "decontamination_controls"),
    value = c(args$seed, args$threads, RNGkind(), paste(predictors[-1], collapse = ","), args$batch_column, args$method,
              "same_inputs_code_environment_seed_threads; no cross-thread-count guarantee", controls))
  conqur <- if (requireNamespace("ConQuR", quietly = TRUE)) utils::packageDescription("ConQuR") else list()
  field <- function(name) if (is.null(conqur[[name]])) "unavailable" else conqur[[name]]
  settings <- rbind(settings, data.frame(key = "selection_distance",
    value = "bray on relative abundance; correction starts only when the baseline batch p < 0.05"))
  if (!is.null(args$depth_summary))
    settings <- rbind(settings, data.frame(key = c("depth_by_batch", "depth_covariate"),
      value = c(args$depth_summary, if (isTRUE(args$depth_covariate)) "added (log_depth = log10 total count)" else "not added")))
  settings <- rbind(settings, data.frame(
    key = c("ConQuR_version", "ConQuR_upstream_sha", "ConQuR_local_patch", "vegan_version", "tuning_permanova"),
    value = c(field("Version"), field("RemoteSha"), field("CMPipelinePatch"),
              as.character(utils::packageVersion("vegan")), "adonis2; bray; by=terms; permutations=999")))
  write.table(settings, file.path(out, "batch_correction_run_settings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  mm <- result$method_matrix
  if (is.null(mm) || (!is.null(result$decision) && grepl("^uncorrected_", result$decision))) {
    methods <- data.frame(prevalence_bin = "all",
                          selection = if (!is.null(result$decision) && grepl("^uncorrected_", result$decision)) result$decision else
                            if (isTRUE(result$success)) "vanilla_requested" else "uncorrected_explicit_fallback",
                          taxa_in_bin = ncol(tab), fit_diagnostics = "not_exposed_by_pinned_dependency")
  } else {
    methods <- as.data.frame(t(mm), stringsAsFactors = FALSE)
    methods$prevalence_bin <- rownames(methods)
    prevalence <- colMeans(tab > 0)
    methods$taxa_in_bin <- vapply(strsplit(methods$prevalence_bin, "-", fixed = TRUE), function(bounds) {
      b <- as.numeric(bounds); sum(prevalence > b[1] & prevalence <= b[2])
    }, integer(1))
    methods$selection <- ifelse(methods$taxa_in_bin == 0, "empty_bin",
                          ifelse(is.na(methods$original), "unknown_selection",
                            ifelse(methods$original == "TRUE", "original_selected", "correction_selected")))
    methods$fit_diagnostics <- "not_exposed_by_pinned_dependency"
    rownames(methods) <- NULL
  }
  write.table(methods, file.path(out, "batch_correction_methods.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("QUALITY: ", q$quality_status, "; biological-signal preservation not assessed; per-taxon fit outcomes unavailable")
}
