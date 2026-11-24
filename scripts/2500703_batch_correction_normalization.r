#!/usr/bin/env Rscript
# ==============================================================================
#  BATCH-CORRECTION & NORMALISATION PIPELINE - Auto-Phased Tune_ConQuR
#  
#  Key improvements:
#  - Phase 1: Fast test (16 combinations, ~10-20 min)
#  - Phase 3: Comprehensive search (72 combinations, ~60-90 min)
#  - PERMANOVA-based stopping criteria (25% R² reduction)
#  - Automatic phase selection based on results






# Example usage (basic):
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --tumor_only

# Example usage (custom columns):
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --batch_column "batch_id" \
#   --covariates "age,gender,BMI" \
#   --type_column "sample_type" \
#   --tumor_value "Tumor,Primary Tumor,Recurrent Tumor" \
#   --tumor_only

# ==============================================================================

# Everything comes from the explicit batch environment lock and its local ConQuR/cqrReg packages.
# The ConQuR upstream commit and local compatibility patch are checked by the module. Nothing is installed at run time: the earlier block installed
# pacman/BiocManager/devtools, an unpinned ConQuR and "DEICODE" (a Python package) on every run,
# which silently changed the software between runs (2026-09-23).
suppressPackageStartupMessages({
  for (pkg in c("tidyverse", "vegan", "ConQuR", "argparse", "zCompositions",
                "doParallel", "foreach", "iterators")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("ERROR: R package '", pkg, "' is missing from the batch-correction environment. ",
           "Rebuild a new prefix from the explicit lock with scripts/setup_batch_env.sh.")
    library(pkg, character.only = TRUE)
  }
})
## RCLR transform only (log ratio to the geometric mean of the non-zero values; the caller sets zero entries
## to 0). This is not robust PCA (RPCA, e.g. DEICODE/Gemelli, which adds matrix completion); PERMANOVA on it
## with Euclidean distance is a sensitivity analysis.
rclr_f <- function(mat){
  t(apply(mat, 1, function(x){
    nz <- x>0; g <- exp(mean(log(x[nz]))); y <- rep(NA, length(x))
    y[nz] <- log(x[nz]/g); y
  }))
}

# ---------- Utility functions -------------------------------------------------
save_table <- function(df, path){
  if (!is.data.frame(df))             
      df <- as.data.frame(df)
  
  # Explicitly preserve column names
  col_names <- colnames(df)
  
  df_out <- df %>%
    rownames_to_column("sample_id")
  
  # Force column names to be written
  write_tsv(df_out, path, progress = FALSE, col_names = TRUE)
  
  cat("  → Saved", nrow(df_out), "rows x", ncol(df_out)-1, "taxa to", basename(path), "\n")
}

# Shared strict readers are sourced from this script's package directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot locate batch script")
source(file.path(dirname(normalizePath(sub("^--file=", "", script_arg))), "batch_contract.R"))
source(file.path(dirname(normalizePath(sub("^--file=", "", script_arg))), "batch_reporting.R"))

# ---------- Select optimal reference batch ------------------------------------
select_optimal_reference_batch <- function(tab, meta, batch_col) {
  
  cat("\n=== SELECTING OPTIMAL REFERENCE BATCH ===\n")
  
  batches <- unique(meta[[batch_col]])
  batch_distances <- matrix(NA, length(batches), length(batches))
  rownames(batch_distances) <- batches
  colnames(batch_distances) <- batches
  
  batch_centroids <- list()
  for (batch in batches) {
    batch_samples <- rownames(meta[meta[[batch_col]] == batch, ])
    batch_data <- tab[batch_samples, , drop = FALSE]
    batch_centroids[[batch]] <- colMeans(batch_data)
  }
  
  for (i in 1:length(batches)) {
    for (j in 1:length(batches)) {
      if (i != j) {
        dist_matrix <- vegdist(rbind(
          batch_centroids[[batches[i]]], 
          batch_centroids[[batches[j]]]
        ), method = "bray")
        batch_distances[i, j] <- as.numeric(dist_matrix)
      } else {
        batch_distances[i, j] <- 0
      }
    }
  }
  
  mean_distances <- rowMeans(batch_distances, na.rm = TRUE)
  
  result_df <- data.frame(
    Batch = names(mean_distances),
    Mean_Distance = mean_distances,
    N_Samples = sapply(batches, function(b) sum(meta[[batch_col]] == b))
  )
  result_df <- result_df[order(result_df$Mean_Distance), ]
  
  cat("\nBatch centrality ranking:\n")
  print(result_df, row.names = FALSE)
  
  top_batches <- result_df$Batch[1:min(3, nrow(result_df))]
  
  cat("\nTop", length(top_batches), "most central batches selected\n")
  
  return(as.character(top_batches))
}

# ---------- PERMANOVA ---------------------------------------------------------
run_permanova <- function(mat, meta, rhs_formula, dist_meth, transformed = FALSE, mode = NULL){
  smp <- rownames(mat)
  dist <- batch_distance(mat, meta, dist_meth, transformed)
  # drop = FALSE: with ONE predictor (batch only -- `--covariates none`, or every covariate
  # filtered out) `[` returned a vector and adonis2 failed "'.' in formula and no 'data' argument".
  pdat <- meta[smp, all.vars(rhs_formula), drop = FALSE]
  # by = "margin" can fail when a strong correction makes the Bray batch sum of squares negative
  # Fall back to sequential terms with batch entered LAST (so its R2 is still the
  # batch effect after covariates), and say so, rather than crashing a successful correction.
  # mode = "terms_batch_last" forces the fallback, so several tables can be compared under one mode.
  evaluation_mode <- "margin"
  terms_last <- function() {
    evaluation_mode <<- "terms_batch_last"
    vars <- all.vars(rhs_formula); vars <- c(setdiff(vars, vars[1]), vars[1])
    adonis2(reformulate(vars, response = "dist"), data = pdat, permutations = 999, by = "terms")
  }
  res <- if (identical(mode, "terms_batch_last")) terms_last() else tryCatch(
    adonis2(dist ~ ., data = pdat, permutations = 999, by = "margin"),
    error = function(e) {
      cat("WARN: PERMANOVA by='margin' failed (", conditionMessage(e), "); using by='terms' with batch last\n")
      terms_last()
    })
  as_tibble(res, rownames = "Term") %>% dplyr::select(Term, Df, R2, `Pr(>F)`) %>%
    mutate(Evaluation = evaluation_mode, Distance = dist_meth, N = nrow(mat),
           Model = paste(deparse(rhs_formula), collapse = " "))
}

# Dispersion check: PERMANOVA R2 can reflect a difference in spread between groups rather
# than in location. betadisper on the same Bray distances, per table and grouping; a failure is recorded, not fatal.
dispersion_check <- function(tables, meta, groups) {
  rows <- list()
  for (tag in names(tables)) for (grp in groups) {
    rows[[length(rows) + 1]] <- tryCatch({
      bd <- betadisper(vegdist(tables[[tag]], "bray"), droplevels(factor(meta[rownames(tables[[tag]]), grp])))
      pt <- permutest(bd, permutations = 999)
      data.frame(Source = tag, Group = grp, F = pt$tab[1, "F"], p = pt$tab[1, "Pr(>F)"], N = nrow(tables[[tag]]),
                 Distance = "bray", Status = "computed")
    }, error = function(e) data.frame(Source = tag, Group = grp, F = NA_real_, p = NA_real_, N = nrow(tables[[tag]]),
                                       Distance = "bray", Status = paste("not computed:", conditionMessage(e))))
  }
  do.call(rbind, rows)
}

# Within-patient permutations: samples from one patient (e.g. TCGA tumour and solid normal) are not independent. R2 does not
# depend on the permutation scheme; for terms that vary within patients the p-value is recomputed with
# permutations restricted within patients. Terms constant within patients keep the free-permutation p only.
patient_restricted_p <- function(mat, meta, rhs_formula, dist_meth, transformed, patient) {
  if (is.null(patient)) return(NULL)
  pt <- patient[rownames(mat)]
  if (anyNA(pt)) return(list(note = "patient missing for some samples; not computed"))
  multi <- sum(table(pt) > 1L)
  if (multi == 0L) return(list(note = "no patient has more than one sample; not computed", n_patients = length(unique(pt)), n_multi = 0L))
  vars <- all.vars(rhs_formula)
  pdat <- meta[rownames(mat), vars, drop = FALSE]
  varies <- vapply(vars, function(v) any(tapply(as.character(pdat[[v]]), pt, function(x) length(unique(x)) > 1L)), logical(1))
  d <- batch_distance(mat, meta, dist_meth, transformed)
  res <- tryCatch(adonis2(d ~ ., data = pdat, by = "margin",
                          permutations = permute::how(blocks = factor(pt), nperm = 999)),
                  error = function(e) NULL)
  p <- setNames(rep(NA_real_, length(vars)), vars)
  if (!is.null(res)) for (v in vars[varies]) if (v %in% rownames(res)) p[[v]] <- res[v, "Pr(>F)"]
  list(p = p, n_patients = length(unique(pt)), n_multi = multi,
       note = if (is.null(res)) "within-patient permutation failed" else "computed")
}

# ---------- Evaluate batch effect (PERMANOVA-based) ---------------------------
evaluate_batch_effect <- function(data, meta, batch_col, covariates, mode = NULL) {
  
  formula <- reformulate(c(batch_col, covariates))
  # Judge the batch effect on relative abundances. Bray on raw counts read a pure
  # sequencing-depth difference between batches as a batch effect (R2 0.15-0.37 with identical composition;
  # synthetic check). A table with a zero-total sample cannot be evaluated.
  data <- as.matrix(data)
  if (any(!is.finite(rowSums(data))) || any(rowSums(data) <= 0))
    return(list(batch_r2 = NA, batch_p = NA, covariate_r2 = NA, mode = NA_character_, full_result = NULL))
  perm_result <- run_permanova(data / rowSums(data), meta, formula, "bray", transformed = TRUE, mode = mode)
  eval_mode <- perm_result$Evaluation[1]
  
  batch_r2 <- perm_result %>% 
    filter(Term == batch_col) %>% 
    pull(R2)
  
  if (length(batch_r2) == 0) {
    warning("Batch term not found in PERMANOVA results")
    return(list(
      batch_r2 = NA,
      batch_p = NA,
      covariate_r2 = NA,
      mode = eval_mode,
      full_result = perm_result
    ))
  }
  
  batch_p <- perm_result %>% 
    filter(Term == batch_col) %>% 
    pull(`Pr(>F)`)
  
  # Covariate R² (for overcorrection check)
  covariate_r2 <- perm_result %>%
    filter(Term %in% covariates) %>%
    pull(R2) %>%
    sum()
  
  return(list(
    batch_r2 = batch_r2,
    batch_p = batch_p,
    covariate_r2 = covariate_r2,
    mode = eval_mode,
    full_result = perm_result
  ))
}

# One PERMANOVA mode for every candidate: candidates are compared only under one PERMANOVA mode. When the evaluations of the
# baseline and the candidates used different modes (margin failed for some), all are re-evaluated with
# sequential terms and batch last. evals and tables are named lists in the same order.
align_eval_modes <- function(evals, tables, meta, batch_col, covariates) {
  modes <- vapply(evals, function(e) if (is.null(e$mode)) NA_character_ else as.character(e$mode), character(1))
  if (length(unique(modes)) <= 1L) return(list(evals = evals, note = ""))
  cat("NOTE: PERMANOVA modes differ (", paste(names(modes), modes, sep = "=", collapse = ", "),
      "); re-evaluating all with terms_batch_last\n")
  evals <- lapply(tables, function(t) evaluate_batch_effect(t, meta, batch_col, covariates, mode = "terms_batch_last"))
  list(evals = evals, note = paste0("; all candidates re-evaluated with terms_batch_last because modes differed (",
                                    paste(names(modes), modes, sep = "=", collapse = ", "), ")"))
}

# ---------- Run single phase --------------------------------------------------
run_phase <- function(phase_num, tab, meta, batch_col, covars, 
                      ref_batch_pool, final_covars) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("PHASE", phase_num, "EXECUTION\n")
  cat(strrep("=", 70), "\n")
  
  # Define parameter pools
  if (phase_num == 1) {
    # PHASE 1: 16 COMBINATIONS (~10-20 minutes)
    batch_pool <- ref_batch_pool[1:min(2, length(ref_batch_pool))]
    logistic_pool <- c(FALSE, TRUE)
    quantile_pool <- c("standard", "lasso")
    lambda_pool <- c(NA, "2p/logn")
    interplt_pool <- c(FALSE)
    
    cat("Strategy: Fast essential coverage\n")
    
  } else {
    # PHASE 3: 72 COMBINATIONS (~60-90 minutes)
    batch_pool <- ref_batch_pool
    logistic_pool <- c(FALSE, TRUE)
    quantile_pool <- c("standard", "lasso")
    lambda_pool <- c(NA, "2p/n", "2p/logn")
    interplt_pool <- c(FALSE, TRUE)
    
    cat("Strategy: Comprehensive grid search\n")
  }
  
  n_combinations <- length(batch_pool) * length(logistic_pool) * 
                    length(quantile_pool) * length(lambda_pool) * 
                    length(interplt_pool)
  
  cat("\nParameter space:\n")
  cat("  Reference batches:", length(batch_pool), "\n")
  cat("  Logistic lasso:   ", length(logistic_pool), "(FALSE, TRUE)\n")
  cat("  Quantile type:    ", length(quantile_pool), "(standard, lasso)\n")
  cat("  Lambda:           ", length(lambda_pool), "\n")
  cat("  Interpolation:    ", length(interplt_pool), "\n")
  cat("\nTotal combinations:", n_combinations, "\n")
  
  if (phase_num == 1) {
    cat("Estimated time:    10-20 minutes\n")
  } else {
    cat("Estimated time:    60-90 minutes\n")
  }
  cat(strrep("=", 70), "\n\n")
  
  start_time <- Sys.time()

  # Run Tune_ConQuR - Keep it simple like the working original
  options(warn = -1)
  registerDoSEQ()

  # CRITICAL: Assign batchid to global environment for foreach workers
  # ConQuR uses foreach internally and workers need access to batchid
  batchid <- factor(meta[[batch_col]])
  assign("batchid", batchid, envir = .GlobalEnv)

  tryCatch({
    # Tune_ConQuR with the single-batch guard on its uncorrected benchmark (batch_contract.R).
    result <- guarded_tune_conqur()(
      tax_tab = tab,
      batchid = batchid,
      covariates = final_covars,
      batch_ref_pool = batch_pool,
      logistic_lasso_pool = logistic_pool,
      quantile_type_pool = quantile_pool,
      simple_match_pool = FALSE,
      lambda_quantile_pool = lambda_pool,
      interplt_pool = interplt_pool,
      frequencyL = 0,
      frequencyU = 1,
      num_core = args$threads
    )
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    options(warn = 0)
    
    cat("\n✓ Phase", phase_num, "completed in", round(elapsed, 1), "minutes\n")
    
    return(list(
      data = as.data.frame(result$tax_final),
      method_matrix = result$method_final,
      elapsed = elapsed,
      phase = phase_num,
      success = TRUE
    ))
    
  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    options(warn = 0)

    cat("\n")
    cat(strrep("!", 70), "\n")
    cat("!!! PHASE", phase_num, "FAILED after", round(elapsed, 1), "minutes !!!\n")
    cat(strrep("!", 70), "\n")
    cat("Error message:\n")
    cat("  ", e$message, "\n")
    cat(strrep("!", 70), "\n\n")

    return(list(
      data = NULL,
      method_matrix = NULL,
      elapsed = elapsed,
      phase = phase_num,
      success = FALSE,
      error = e$message
    ))
  })
}

# ---------- Vanilla ConQuR (default, non-tuned) --------------------------------
run_vanilla_conqur <- function(tab_raw, meta, batch_col, covars, corr_dir) {

  cat("\n", strrep("=", 70), "\n")
  cat("VANILLA ConQuR MODE (non-tuned)\n")
  cat(strrep("=", 70), "\n\n")

  # The cohort and requested covariates were fixed once before dispatch.
  tab <- tab_raw
  valid_covars <- covars
  final_covars <- meta[, covars, drop = FALSE]

  # Baseline evaluation
  baseline <- evaluate_batch_effect(tab, meta, batch_col, valid_covars)
  assign("decision_baseline", baseline, envir = .GlobalEnv)   # values the correction decision used, for the report
  cat(sprintf("Baseline Batch R2: %.4f (p = %.4f)\n", baseline$batch_r2, baseline$batch_p))
  d0 <- decide_correction(baseline$batch_r2, numeric(0), baseline_p = baseline$batch_p)
  if (d0$decision %in% c("uncorrected_not_evaluable", "uncorrected_no_batch_effect"))
    return(uncorrected_by_decision(tab, baseline$batch_r2, d0))

  # Select reference batch
  ref_batch_pool <- select_optimal_reference_batch(tab, meta, batch_col)
  batch_ref <- ref_batch_pool[1]
  cat("Using reference batch:", batch_ref, "\n\n")

  # Run vanilla ConQuR
  batchid <- factor(meta[[batch_col]])
  assign("batchid", batchid, envir = .GlobalEnv)

  start_time <- Sys.time()
  options(warn = -1)
  registerDoSEQ()

  tryCatch({
    result <- ConQuR(
      tax_tab = tab,
      batchid = batchid,
      covariates = final_covars,
      batch_ref = batch_ref,
      logistic_lasso = FALSE,
      quantile_type = "lasso",
      interplt = FALSE,
      num_core = args$threads
    )

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    options(warn = 0)

    corrected <- as.data.frame(result)

    # Evaluate
    post_eval <- evaluate_batch_effect(corrected, meta, batch_col, valid_covars)
    delta_r2 <- baseline$batch_r2 - post_eval$batch_r2
    r2_reduction <- delta_r2 / baseline$batch_r2

    cat(sprintf("\nVanilla ConQuR completed in %.1f minutes\n", elapsed))
    cat(sprintf("Before: Batch R2 = %.4f\n", baseline$batch_r2))
    cat(sprintf("After:  Batch R2 = %.4f\n", post_eval$batch_r2))
    cat(sprintf("Reduction: %.1f%%\n", report_reduction_percent(baseline$batch_r2, baseline$batch_r2 - delta_r2)))
    d <- decide_correction(baseline$batch_r2, c(vanilla = post_eval$batch_r2))
    if (d$decision != "corrected") return(uncorrected_by_decision(tab, baseline$batch_r2, d, elapsed))

    return(list(
      decision = d$decision, decision_reason = d$reason,
      data = corrected,
      method_matrix = NULL,
      phase_used = "vanilla",
      baseline_r2 = baseline$batch_r2,
      final_r2 = post_eval$batch_r2,
      r2_reduction = r2_reduction,
      elapsed = elapsed,
      success = TRUE
    ))

  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    options(warn = 0)
    cat("\nVanilla ConQuR FAILED after", round(elapsed, 1), "minutes\n")
    cat("Error:", e$message, "\n")
    return(list(data = tab, phase_used = "vanilla", success = FALSE, message = e$message))
  })
}

# ---------- Decision outcome that keeps the uncorrected table -------------------
# A deliberate outcome, not a failure: the run continues, the table is
# published as UNCORRECTED_by_decision.tsv (never as ConQuR_tuned.tsv) and the reason is recorded.
uncorrected_by_decision <- function(tab, baseline_r2, d, elapsed = 0) {
  cat("✓ DECISION: keep the uncorrected table (", d$decision, "): ", d$reason, "\n", sep = "")
  list(data = tab, method_matrix = NULL, phase_used = 0, baseline_r2 = baseline_r2,
       final_r2 = baseline_r2, r2_reduction = NA_real_, elapsed = elapsed, success = TRUE,
       decision = d$decision, decision_reason = d$reason)
}

# ---------- Auto-phase workflow -----------------------------------------------
run_auto_phase <- function(tab_raw, meta, batch_col, covars,
                            corr_dir, r2_threshold = 0.25) {

  cat("\n", strrep("=", 70), "\n")
  cat("AUTO-PHASE MODE: PERMANOVA-BASED STRATEGY\n")
  cat(strrep("=", 70), "\n")
  cat("Decision rule:\n")
  cat("  IF R² reduction ≥", r2_threshold * 100, "% → STOP at Phase 1\n", sep="")
  cat("  ELSE → Run Phase 3 and compare results\n")
  cat(strrep("=", 70), "\n\n")

  # Identical cohort/model for correction and every raw/corrected comparison.
  tab <- tab_raw
  if (nrow(meta) < 20) {
    return(list(data = tab, phase_used = 0, success = FALSE,
                message = "Insufficient samples: tuned mode requires at least 20 eligible samples"))
  }
  valid_covars <- covars
  final_covars <- meta[, covars, drop = FALSE]

  # Baseline evaluation
  cat("\n", strrep("-", 70), "\n")
  cat("BASELINE EVALUATION (Before correction)\n")
  cat(strrep("-", 70), "\n")

  baseline <- evaluate_batch_effect(tab, meta, batch_col, valid_covars)
  assign("decision_baseline", baseline, envir = .GlobalEnv)   # values the correction decision used, for the report
  
  cat(sprintf("Batch R²:     %.4f (p = %.4f)\n", 
              baseline$batch_r2, baseline$batch_p))
  if (!is.na(baseline$covariate_r2) && baseline$covariate_r2 > 0) {
    cat(sprintf("Covariate R²: %.4f\n", baseline$covariate_r2))
  }
  cat(strrep("-", 70), "\n\n")
  d0 <- decide_correction(baseline$batch_r2, numeric(0), baseline_p = baseline$batch_p)
  if (d0$decision %in% c("uncorrected_not_evaluable", "uncorrected_no_batch_effect"))
    return(uncorrected_by_decision(tab, baseline$batch_r2, d0))
  
  # Select reference batches
  ref_batch_pool <- select_optimal_reference_batch(tab, meta, batch_col)
  
  # === PHASE 1 ===
  cat("\n>>> STARTING PHASE 1 <<<\n")
  phase1_result <- run_phase(1, tab, meta, batch_col, valid_covars, 
                              ref_batch_pool, final_covars)
  
  if (!phase1_result$success) {
    cat("\n")
    cat(strrep("!", 70), "\n")
    cat("!!! ConQuR PHASE 1 FAILED - RETURNING UNCORRECTED DATA !!!\n")
    cat(strrep("!", 70), "\n\n")
    return(list(
      data = tab,
      phase_used = 0,
      success = FALSE,
      message = phase1_result$error
    ))
  }
  
  # Evaluate Phase 1
  cat("\n>>> EVALUATING PHASE 1 <<<\n")
  phase1_eval <- evaluate_batch_effect(phase1_result$data, meta, 
                                        batch_col, valid_covars)
  al <- align_eval_modes(list(baseline = baseline, phase1 = phase1_eval),
                         list(baseline = tab, phase1 = phase1_result$data), meta, batch_col, valid_covars)
  baseline <- al$evals$baseline; phase1_eval <- al$evals$phase1; mode_note <- al$note
  
  delta_r2 <- baseline$batch_r2 - phase1_eval$batch_r2
  r2_reduction <- delta_r2 / baseline$batch_r2
  
  cat("\n", strrep("-", 70), "\n")
  cat("PHASE 1 EVALUATION\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("Before:       Batch R² = %.4f (p = %.4f)\n", 
              baseline$batch_r2, baseline$batch_p))
  cat(sprintf("After:        Batch R² = %.4f (p = %.4f)\n", 
              phase1_eval$batch_r2, phase1_eval$batch_p))
  cat(sprintf("Change:       ΔR² = %.4f\n", delta_r2))
  cat(sprintf("Reduction:    %.1f%%\n", report_reduction_percent(baseline$batch_r2, baseline$batch_r2 - delta_r2)))
  
  # Overcorrection check
  if (!is.na(baseline$covariate_r2) && baseline$covariate_r2 > 0) {
    covar_reduction <- (baseline$covariate_r2 - phase1_eval$covariate_r2) / 
                       baseline$covariate_r2
    if (isTRUE(covar_reduction > 0.2)) {
      cat(sprintf("⚠ WARNING:    Covariate R² reduced by %.1f%% (possible overcorrection)\n", 
                  covar_reduction * 100))
    }
  }
  cat(strrep("-", 70), "\n\n")
  
  # Decision point
  if (isTRUE(is.finite(phase1_eval$batch_r2) && phase1_eval$batch_r2 < baseline$batch_r2 && r2_reduction >= r2_threshold)) {
    cat("✓ DECISION: STOP at Phase 1\n")
    cat(sprintf("  R² reduction (%.1f%%) meets threshold (≥%.1f%%)\n", 
                report_reduction_percent(baseline$batch_r2, baseline$batch_r2 - delta_r2), r2_threshold * 100))
    cat("  Phase 3 not needed.\n\n")
    
    # Save Phase 1 result
    final_result <- list(
      data = phase1_result$data,
      method_matrix = phase1_result$method_matrix,
      phase_used = 1,
      baseline_r2 = baseline$batch_r2,
      final_r2 = phase1_eval$batch_r2,
      r2_reduction = r2_reduction,
      elapsed = phase1_result$elapsed,
      success = TRUE
    )
    d <- decide_correction(baseline$batch_r2, c(phase1 = phase1_eval$batch_r2), r2_threshold)
    final_result$decision <- d$decision; final_result$decision_reason <- paste0(d$reason, mode_note)
    
  } else {
    # Need Phase 3
    if (!isTRUE(phase1_eval$batch_r2 < baseline$batch_r2)) {
      cat("⚠ DECISION: RUN PHASE 3\n")
      cat("  Phase 1 did not improve batch effect (R² unchanged or worse)\n")
    } else {
      cat("→ DECISION: RUN PHASE 3\n")
      cat(sprintf("  Phase 1 R² reduction (%.1f%%) below threshold (%.1f%%)\n", 
                  report_reduction_percent(baseline$batch_r2, baseline$batch_r2 - delta_r2), r2_threshold * 100))
      cat("  Phase 3 may find better parameters.\n")
    }
    
    cat("\nWaiting 5 seconds before starting Phase 3...\n")
    Sys.sleep(5)
    
    # === PHASE 3 ===
    cat("\n>>> STARTING PHASE 3 <<<\n")
    phase3_result <- run_phase(3, tab, meta, batch_col, valid_covars, 
                                ref_batch_pool, final_covars)
    
    if (!phase3_result$success) {
      d <- decide_correction(baseline$batch_r2, c(phase1 = phase1_eval$batch_r2), r2_threshold)
      if (d$decision != "corrected") {
        final_result <- uncorrected_by_decision(tab, baseline$batch_r2, d, phase1_result$elapsed)
        final_result$note <- "Phase 3 failed; Phase 1 did not improve"
      } else {
      cat("\nPhase 3 failed. Using Phase 1 result.\n")
      final_result <- list(
        decision = d$decision, decision_reason = paste0(d$reason, mode_note),
        data = phase1_result$data,
        method_matrix = phase1_result$method_matrix,
        phase_used = 1,
        baseline_r2 = baseline$batch_r2,
        final_r2 = phase1_eval$batch_r2,
        r2_reduction = r2_reduction,
        elapsed = phase1_result$elapsed,
        success = TRUE,
        note = "Phase 3 failed, used Phase 1"
      )
      }
    } else {
      # Evaluate Phase 3
      cat("\n>>> EVALUATING PHASE 3 <<<\n")
      phase3_eval <- evaluate_batch_effect(phase3_result$data, meta, 
                                            batch_col, valid_covars)
      al <- align_eval_modes(list(baseline = baseline, phase1 = phase1_eval, phase3 = phase3_eval),
                             list(baseline = tab, phase1 = phase1_result$data, phase3 = phase3_result$data),
                             meta, batch_col, valid_covars)
      baseline <- al$evals$baseline; phase1_eval <- al$evals$phase1; phase3_eval <- al$evals$phase3
      if (nzchar(al$note)) mode_note <- al$note
      delta_r2 <- baseline$batch_r2 - phase1_eval$batch_r2
      
      phase3_delta <- baseline$batch_r2 - phase3_eval$batch_r2
      phase3_reduction <- phase3_delta / baseline$batch_r2
      
      cat("\n", strrep("-", 70), "\n")
      cat("PHASE 3 EVALUATION\n")
      cat(strrep("-", 70), "\n")
      cat(sprintf("Baseline:     Batch R² = %.4f\n", baseline$batch_r2))
      cat(sprintf("Phase 1:      Batch R² = %.4f (%.1f%% reduction)\n", 
                  phase1_eval$batch_r2, report_reduction_percent(baseline$batch_r2, baseline$batch_r2 - delta_r2)))
      cat(sprintf("Phase 3:      Batch R² = %.4f (%.1f%% reduction)\n", 
                  phase3_eval$batch_r2, report_reduction_percent(baseline$batch_r2, phase3_eval$batch_r2)))
      cat(strrep("-", 70), "\n\n")
      
      # Compare and select: neither NA nor a candidate worse than the baseline can be chosen.
      d <- decide_correction(baseline$batch_r2, c(phase1 = phase1_eval$batch_r2, phase3 = phase3_eval$batch_r2), r2_threshold)
      d$reason <- paste0(d$reason, mode_note)
      if (d$decision != "corrected") {
        final_result <- uncorrected_by_decision(tab, baseline$batch_r2, d, phase1_result$elapsed + phase3_result$elapsed)
      } else if (identical(d$choice, "phase3")) {
        improvement <- ((phase1_eval$batch_r2 - phase3_eval$batch_r2) / 
                       phase1_eval$batch_r2) * 100
        cat("✓ DECISION: Use Phase 3 result\n")
        cat(sprintf("  Phase 3 achieved %.1f%% lower Batch R² than Phase 1\n", 
                    improvement))
        
        final_result <- list(
          decision = d$decision, decision_reason = d$reason,
          data = phase3_result$data,
          method_matrix = phase3_result$method_matrix,
          phase_used = 3,
          baseline_r2 = baseline$batch_r2,
          final_r2 = phase3_eval$batch_r2,
          r2_reduction = phase3_reduction,
          elapsed = phase1_result$elapsed + phase3_result$elapsed,
          success = TRUE,
          phase1_r2 = phase1_eval$batch_r2
        )
        
      } else {
        cat("✓ DECISION: Use Phase 1 result\n")
        cat("  Phase 3 did not improve upon Phase 1\n")
        
        final_result <- list(
          decision = d$decision, decision_reason = d$reason,
          data = phase1_result$data,
          method_matrix = phase1_result$method_matrix,
          phase_used = 1,
          baseline_r2 = baseline$batch_r2,
          final_r2 = phase1_eval$batch_r2,
          r2_reduction = r2_reduction,
          elapsed = phase1_result$elapsed + phase3_result$elapsed,
          success = TRUE,
          note = "Phase 1 was better",
          phase3_r2 = phase3_eval$batch_r2
        )
      }
    }
  }
  
  cat("\n", strrep("=", 70), "\n")
  cat("AUTO-PHASE COMPLETED\n")
  cat(strrep("=", 70), "\n")
  cat("Phase used:       ", final_result$phase_used, "\n", sep="")
  cat(sprintf("Baseline R²:      %.4f\n", final_result$baseline_r2))
  cat(sprintf("Final R²:         %.4f\n", final_result$final_r2))
  cat(sprintf("Total reduction:  %.1f%%\n", report_reduction_percent(final_result$baseline_r2, final_result$final_r2)))
  cat(sprintf("Total time:       %.1f minutes\n", final_result$elapsed))
  cat(strrep("=", 70), "\n\n")

  return(final_result)
}

# ---------- PCoA plots --------------------------------------------------------
generate_conqur_plots <- function(tab, meta, batch_col, title, out_pref){
  smp <- intersect(rownames(tab), rownames(meta))
  fac <- meta[smp, batch_col]
  png(paste0(out_pref,"_bray.png"), 8,7,'in',res=300)
  Plot_PCoA(TAX = tab[smp, ], factor = fac,
            main = paste(title,", Bray-Curtis"))
  dev.off()
  png(paste0(out_pref,"_aitch.png"), 8,7,'in',res=300)
  Plot_PCoA(TAX = tab[smp, ], factor = fac, dissimilarity = "Aitch",
            main = paste(title,", Aitchison"))
  dev.off()
}

# ---------- Normalize ---------------------------------------------------------
normalize_all <- function(mat, stub){
  matrix_ok(mat, counts = TRUE)
  if (any(rowSums(mat) == 0)) stop("Cannot normalize zero-total analysis samples")
  dir.create(dirname(stub), showWarnings = FALSE, recursive = TRUE)
  paths <- c(CLR_pseudo = paste0(stub, "_clr_pseudo.tsv"),
             CLR_CZM = paste0(stub, "_clr_czm.tsv"), RCLR = paste0(stub, "_rclr.tsv"))
  if (any(file.exists(paths))) stop("Existing normalization outputs: use a fresh output directory")
  x <- as.matrix(mat)
  out <- list()
  for (method in names(paths)) {
    out[[method]] <- tryCatch({
      result <- switch(method,
        CLR_pseudo = {
          pseudo <- max(as.numeric(quantile(x[x > 0], 0.05)) * 0.5, 0.1)
          decostand(x + pseudo, "clr")
        },
        CLR_CZM = {
          # Strictly positive input needs no zero replacement.
          imputed <- if (any(x == 0)) cmultRepl(x, method = "CZM", label = 0,
                        output = "p-counts", z.warning = 1, z.delete = FALSE,
                        suppress.print = TRUE) else x
          decostand(imputed, "clr")
        },
        RCLR = {
          y <- rclr_f(x)
          # Algorithmic zeros only; do not conceal unrelated NA/NaN failures.
          y[x == 0] <- 0
          dimnames(y) <- dimnames(x)
          y
        })
      result <- as.data.frame(result, check.names = FALSE)
      matrix_ok(result)
      assert_same_matrix_ids(result, mat)
      result
    }, error = function(e) {
      reason <- paste0("Required normalization ", method, " failed: ", conditionMessage(e))
      write_status("failed", reason)
      stop(reason, call. = FALSE)
    })
  }
  # Validate the complete requested set before writing any method for this source.
  for (method in names(paths)) save_table(out[[method]], paths[[method]])
  out
}

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

parser <- ArgumentParser()
parser$add_argument("--otu",    required = TRUE)
parser$add_argument("--meta",   required = TRUE)
parser$add_argument("--prefix", required = TRUE)
parser$add_argument("--tumor_only", action = "store_true")
parser$add_argument("--seed", type = "integer", default = 1)
parser$add_argument("--threads", type = "integer", default = 2)
parser$add_argument("--phase", default = "auto",
                   help = "Phase mode: 'auto', '1', or '3'")
parser$add_argument("--r2_threshold", type = "double", default = 0.25,
                   help = "R2 reduction threshold for stopping at Phase 1 (default: 0.25)")
parser$add_argument("--batch_column", type = "character", default = "shipment_batch",
                   help = "Name of the column in metadata that specifies batch (default: shipment_batch)")
parser$add_argument("--covariates", type = "character", default = "tcga_patient_age",
                   help = "Comma-separated list of covariate column names (default: tcga_patient_age)")
parser$add_argument("--type_column", type = "character", default = "Type",
                   help = "Name of the column in metadata that specifies sample type for tumor filtering (default: Type)")
parser$add_argument("--tumor_value", type = "character", default = "Tumor",
                   help = "Comma-separated list of values to recognize as tumor samples (default: Tumor). Case-insensitive matching.")
parser$add_argument("--taxon_column", type = "character", default = "clade_name",
                    help = "OTU-table column holding taxon names (the pipeline passes params.taxon_column)")
parser$add_argument("--sample_id_column", type = "character", default = "sample_id",
                    help = "metadata column holding sample IDs (the pipeline passes decontam's)")
parser$add_argument("--allow_uncorrected", action = "store_true", default = FALSE,
                    help = "on failure, continue with raw counts saved as UNCORRECTED_passthrough.tsv")
parser$add_argument("--control_values", type = "character", default = NULL,
                    help = "comma-separated decontamination control labels in --type_column; excluded from correction unless --include_controls")
parser$add_argument("--patient_column", type = "character", default = NULL,
                    help = "metadata column naming each sample's patient; restricts PERMANOVA permutations within patients")
parser$add_argument("--effect_column", type = "character", default = NULL,
                    help = "metadata column with two values; effect sizes before/after correction are also reported for it")
parser$add_argument("--patient_from_tcga_barcode", action = "store_true", default = FALSE,
                    help = "use the first 12 characters of the sample ID (TCGA-XX-XXXX) as the patient")
parser$add_argument("--no_protect_sample_type", action = "store_true", default = FALSE,
                    help = "do not add sample type as a protected covariate when several types remain (recorded)")
parser$add_argument("--depth_covariate", action = "store_true", default = FALSE,
                    help = "add log10 total count per sample (log_depth) as a covariate (off by default)")
parser$add_argument("--include_controls", action = "store_true", default = FALSE,
                    help = "keep decontamination controls in the correction cohort (a declared analysis choice)")
parser$add_argument("--method", type = "character", default = "tune",
                   help = "Correction method: 'tune' (default, auto-phased Tune_ConQuR) or 'vanilla' (default ConQuR)")
args <- parser$parse_args()
configure_batch_rng(args$seed, args$threads)
cat("Reproducibility: seed=", args$seed, "; threads=", args$threads, "; RNG=", paste(RNGkind(), collapse=","), "\n", sep="")
cat("R2 reduction NA% means non-interpretable; see batch_correction_quality.tsv. Quality is descriptive, not an efficacy pass.\n")
if (!(args$method %in% c("tune", "vanilla")))
  stop("ERROR: --method must be 'tune' or 'vanilla' (got '", args$method, "')")

reject_existing_batch_outputs(dirname(args$prefix))
dir.create(dirname(args$prefix), showWarnings = FALSE, recursive = TRUE)
corr_dir <- file.path(dirname(args$prefix), "corrected"); dir.create(corr_dir, showWarnings=FALSE)
# One status file per run, whatever happens: corrected | failed, with the reason and sample count.
write_status <- function(status, reason) {
  n <- if (exists("otu_raw", inherits = TRUE)) nrow(get("otu_raw", inherits = TRUE)) else NA
  n_dropped <- if (exists("dropped_empty_samples", inherits = TRUE)) length(get("dropped_empty_samples", inherits = TRUE)) else NA
  write.table(data.frame(status = status, reason = reason, n_samples = n, n_empty_samples_dropped = n_dropped),
              file.path(corr_dir, "BATCH_CORRECTION_STATUS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
norm_dir <- file.path(dirname(args$prefix), "normalized"); dir.create(norm_dir, showWarnings=FALSE)
plot_dir <- file.path(dirname(args$prefix), "pcoa_plots"); dir.create(plot_dir, showWarnings=FALSE)

# Parse covariates (comma-separated string to vector)
covariate_vector <- trimws(unlist(strsplit(args$covariates, ",")))
# "none" (or empty) stops the run: ConQuR needs >= 1 covariate. A named covariate
# missing from the metadata also stops the run below. A covariate that is constant within any
# batch also stops the run (prepare_batch_cohort in batch_contract.R): the requested model is
# never changed silently.
if (length(covariate_vector) == 0 || identical(tolower(covariate_vector), "none")) {
  # ConQuR builds its model from covariates + batch and cannot run on batch alone.
  # Decided 2026-09-23 (Ammal): fail loudly rather than guess or pass raw counts off as corrected.
  write_status("failed", "no covariates requested (--covariates none); ConQuR needs >= 1 covariate that varies within batches")
  stop("ERROR: ConQuR needs >= 1 covariate that varies within batches; --covariates none cannot be batch-corrected.")
}
if (!identical(args$phase, "auto")) {
  # run_auto_phase always runs phase 1 then (if needed) phase 3; forcing a phase was never wired.
  stop("ERROR: --phase '", args$phase, "' is not implemented; only 'auto' is.")
}
cat("INFO: Using batch column:", args$batch_column, "\n")
cat("INFO: Using covariates:", paste(covariate_vector, collapse = ", "), "\n")
if (args$tumor_only) {
  # Parse tumor values (comma-separated string to vector)
  tumor_value_vector <- trimws(unlist(strsplit(args$tumor_value, ",")))
  cat("INFO: Using type column for filtering:", args$type_column, "\n")
  cat("INFO: Tumor values to match:", paste(tumor_value_vector, collapse = ", "), "\n")
}

# Validate all IDs first, then select one cohort for every downstream operation.
loaded <- load_and_prepare_data(args$otu, args$meta)
data <- prepare_batch_cohort(loaded, args, covariate_vector)
otu_raw <- data$otu
meta <- data$meta
covariate_vector <- data$covariates   # may now include sample_type_protected and log_depth
args$depth_summary <- data$depth_summary
# Patient of each sample, from --patient_column or, with --patient_from_tcga_barcode, the first 12
# characters of the sample ID (TCGA-XX-XXXX). NULL when neither is given.
patient_vector <- if (!is.null(args$patient_column) && nzchar(args$patient_column)) {
  if (!args$patient_column %in% names(loaded$meta)) stop("Declared patient column not found: ", args$patient_column)
  setNames(trimws(as.character(loaded$meta[rownames(meta), args$patient_column])), rownames(meta))
} else if (isTRUE(args$patient_from_tcga_barcode)) setNames(substr(rownames(meta), 1, 12), rownames(meta)) else NULL

predictors  <- c(args$batch_column, covariate_vector)
rhs_formula <- reformulate(predictors)

## Run ConQuR with auto-phase
tables <- list(raw = otu_raw)

cat("\n", strrep("=", 70), "\n")
cat("STARTING BATCH CORRECTION WITH AUTO-PHASE STRATEGY\n")
cat(strrep("=", 70), "\n")
cat("Mode:", args$phase, "\n")
cat("R² threshold:", args$r2_threshold * 100, "%\n")
cat(strrep("=", 70), "\n\n")

# Check cache
cache_file <- file.path(corr_dir, "ConQuR_tuned.tsv")
param_file <- file.path(corr_dir, "ConQuR_tuned_parameters.txt")

{
  # Dispatch based on --method argument
  if (args$method == "vanilla") {
    cat("INFO: Using VANILLA ConQuR method\n")
    result <- run_vanilla_conqur(
      otu_raw, meta, args$batch_column,
      predictors[-1],
      corr_dir = corr_dir
    )
  } else {
    # Default: auto-phased Tune_ConQuR
    result <- run_auto_phase(
      otu_raw, meta, args$batch_column,
      predictors[-1],
      corr_dir = corr_dir,
      r2_threshold = args$r2_threshold
    )
  }
  
  result$sample_type_note <- data$sample_type_note
  result$confounding_note <- data$confounding_note
  decided_uncorrected <- isTRUE(result$success) && !is.null(result$decision) && grepl("^uncorrected_", result$decision)
  if (decided_uncorrected) {
    # Deliberate outcome, not a failure: continue with the input table under an honest name.
    tables$UNCORRECTED_by_decision <- otu_raw
    save_table(otu_raw, file.path(corr_dir, "UNCORRECTED_by_decision.tsv"))
    write_status("in_progress", paste0(result$decision, ": ", result$decision_reason))
    sink(param_file)
    cat("Auto-Phase Tune_ConQuR Results\n")
    cat("==============================\n\n")
    cat("Decision: ", result$decision, "\n", sep = "")
    cat("Reason:   ", result$decision_reason, "\n", sep = "")
    cat("\nThe uncorrected table is used downstream (UNCORRECTED_by_decision.tsv).\n")
    sink()
  } else if (result$success) {
    matrix_ok(result$data, counts = TRUE)
    assert_same_matrix_ids(result$data, otu_raw)
    tables$ConQuR_tuned <- result$data
    write_status("in_progress", "ConQuR returned; required downstream outputs not yet complete")
    
    # Save results
    save_table(result$data, cache_file)
    
    # Save parameters
    sink(param_file)
    cat("Auto-Phase Tune_ConQuR Results\n")
    cat("==============================\n\n")
    cat("Phase used:       ", result$phase_used, "\n", sep="")
    cat("Total time:       ", round(result$elapsed, 1), " minutes\n", sep="")
    cat("\nQuality Metrics:\n")
    cat(sprintf("  Baseline Batch R²:     %.4f\n", result$baseline_r2))
    cat(sprintf("  Final Batch R²:        %.4f\n", result$final_r2))
    cat(sprintf("  R² reduction:          %.1f%%\n", report_reduction_percent(result$baseline_r2, result$final_r2)))
    if (!is.null(result$phase1_r2)) {
      cat(sprintf("  Phase 1 Batch R²:      %.4f\n", result$phase1_r2))
    }
    if (!is.null(result$phase3_r2)) {
      cat(sprintf("  Phase 3 Batch R²:      %.4f\n", result$phase3_r2))
    }
    cat("\nSelected parameters by taxa frequency:\n")
    cat("=====================================\n")
    print(result$method_matrix)
    cat("\n")
    if (!is.null(result$note)) {
      cat("Note:", result$note, "\n")
    }
    sink()
    
    cat("\n✓ Parameters saved to:", param_file, "\n")
    
  } else {
    # ==================== BATCH CORRECTION FAILED ====================
    cat("\n")
    cat(strrep("!", 70), "\n")
    cat("!!! ERROR: BATCH CORRECTION (ConQuR) FAILED !!!\n")
    cat(strrep("!", 70), "\n")
    cat("Reason:", ifelse(is.null(result$message), "Unknown error", result$message), "\n")
    cat(strrep("!", 70), "\n\n")

    # Never publish raw counts as ConQuR_tuned.tsv: record the failure and stop,
    # unless --allow_uncorrected, in which case the raw table is saved under an honest name.
    reason <- ifelse(is.null(result$message), "Unknown error", result$message)
    write_status("failed", reason)
    if (!isTRUE(args$allow_uncorrected)) {
      stop("ERROR: batch correction failed: ", reason,
           ". No corrected table was written (see corrected/BATCH_CORRECTION_STATUS.tsv). ",
           "Pass --allow_uncorrected to continue with UNCORRECTED data under that name.")
    }
    tables$UNCORRECTED_passthrough <- otu_raw
    save_table(otu_raw, file.path(corr_dir, "UNCORRECTED_passthrough.tsv"))

    # Save failure info to parameters file
    sink(param_file)
    cat("Auto-Phase Tune_ConQuR Results\n")
    cat("==============================\n\n")
    cat("Status: FAILED\n")
    cat("Reason:", ifelse(is.null(result$message), "Unknown error", result$message), "\n")
    cat("\nUsing uncorrected data for downstream analysis.\n")
    sink()

    cat("WARNING: --allow_uncorrected: saved UNCORRECTED data as UNCORRECTED_passthrough.tsv\n\n")
  }
}

cat("\n", strrep("=", 70), "\n")
cat("BATCH CORRECTION COMPLETED\n")
cat(strrep("=", 70), "\n\n")

## PERMANOVA summary
sum_path <- file.path(dirname(args$prefix), "permanova_summary.tsv")
if (file.exists(sum_path)) stop("Existing PERMANOVA output: use a fresh output directory")

## Main loop
permsum <- list()
for (tag in names(tables)){
  mat <- tables[[tag]]

  # A plot is supporting output: if PCoA cannot be drawn (e.g. no positive eigenvalue when all
  # distances are zero), record why and keep the decision and evaluation records (2026-11-17).
  plot_error <- tryCatch({
    generate_conqur_plots(mat, meta, args$batch_column, tag, file.path(plot_dir, tag)); NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(plot_error)) {
    writeLines(paste("PCoA plot skipped:", plot_error), file.path(plot_dir, paste0(tag, "_PLOT_SKIPPED.txt")))
    cat("WARN: PCoA plot skipped for", tag, ":", plot_error, "\n")
  }

  pp <- patient_restricted_p(mat, meta, rhs_formula, "bray", FALSE, patient_vector)
  permsum[[length(permsum)+1]] <- run_permanova(
        mat, meta, rhs_formula, "bray") %>%
        mutate(Source = tag, Transform = "counts",
               p_within_patient = if (!is.null(pp$p)) unname(pp$p[Term]) else NA_real_,
               Patients = if (is.null(pp$n_patients)) NA_integer_ else pp$n_patients,
               Patients_multi = if (is.null(pp$n_multi)) NA_integer_ else pp$n_multi,
               Patient_note = if (is.null(pp)) "no patient column declared" else pp$note)

  norm_res <- normalize_all(mat,
        file.path(norm_dir, paste0(basename(args$prefix),"_",tag)))

  if (length(norm_res) > 0) {
    for (tr in names(norm_res)){
      if (!is.null(norm_res[[tr]])) {
        permsum[[length(permsum)+1]] <- run_permanova(
              norm_res[[tr]], meta, rhs_formula, "euclidean", transformed = TRUE) %>%
              mutate(Source = tag, Transform = tr)
      }
    }
  } else {
    cat("WARN: All normalization failed for", tag, "\n")
  }
}

bind_rows(permsum) %>% save_table(sum_path)
# Dispersion check: PERMANOVA R2 can reflect a difference in spread between groups rather
# than in location. betadisper on the same Bray distances, per table, for the batch and (when protected) the
# sample type; a failure is recorded, not fatal.
write.table(dispersion_check(tables, meta, intersect(c(args$batch_column, "sample_type_protected"), names(meta))),
            file.path(dirname(args$prefix), "batch_correction_dispersion.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write_effect_sizes(tables, meta, args$type_column, trimws(unlist(strsplit(args$tumor_value, ","))),
                   file.path(dirname(args$prefix), "batch_correction_effect_sizes.tsv"))
if (!is.null(args$effect_column) && nzchar(args$effect_column))
  write_named_effect_sizes(tables, meta, args$effect_column, covariate_vector, args$type_column,
                           file.path(dirname(args$prefix), "batch_correction_effect_sizes.tsv"))
write_batch_records(result, otu_raw, predictors, args)
if ("ConQuR_tuned" %in% names(tables))
  write_status("corrected", "ConQuR and all required normalization/evaluation outputs completed")
if ("UNCORRECTED_by_decision" %in% names(tables))
  write_status(result$decision, result$decision_reason)

cat("\n", strrep("=", 70), "\n")
cat("PIPELINE COMPLETED\n")
cat(strrep("=", 70), "\n")
cat("Corrected data:   ", cache_file, "\n", sep="")
cat("Parameters:       ", param_file, "\n", sep="")
cat("PERMANOVA results:", sum_path, "\n")
cat("PCoA plots:       ", plot_dir, "\n", sep="")
cat("\nPlease review PERMANOVA R² values to assess batch correction quality.\n")
cat(strrep("=", 70), "\n\n")

message("DONE")

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# Auto mode (basic, recommended)
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --tumor_only

# Custom columns and covariates
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --batch_column "batch_id" \
#   --covariates "age,gender,BMI" \
#   --type_column "sample_type" \
#   --tumor_value "Tumor,Primary Tumor,Recurrent Tumor" \
#   --tumor_only

# Custom R² threshold (default 25%)
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --r2_threshold 0.30  # 30% reduction required

# (--phase accepts only 'auto': manual phases 1/3 were removed; the pipeline refuses others at launch)

# Multiple tumor types (case-insensitive)
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --tumor_value "Tumor,Primary Tumor,Recurrent Tumor,Metastatic" \
#   --tumor_only

# ==============================================================================
