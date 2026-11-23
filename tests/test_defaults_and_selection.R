#!/usr/bin/env Rscript
# Defaults and selection: decontamination controls leave the correction cohort
# by default; mixed sample types without sample type among the covariates are recorded, not changed
# or stopped; automatic selection never stops or selects on 0/NA and keeps the uncorrected table
# when no evaluable candidate improves on the baseline.
# Needs only base R plus readr/dplyr/tibble. Fitting and PERMANOVA are injected, not run.
# Usage: Rscript tests/test_defaults_and_selection.R <repo_root> <scratch_dir>
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})
argv <- commandArgs(TRUE); root <- argv[1]; out <- argv[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "scripts/batch_contract.R"))
source(file.path(root, "scripts/batch_reporting.R"))
code <- parse(file.path(root, "scripts/2500703_batch_correction_normalization.r"))
for (expr in code) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) eval(expr)
}
passed <- 0L
test <- function(name, expr) { force(expr); passed <<- passed + 1L; cat("PASS", name, "\n") }

# ---------------------------------------------------------------- decide_correction
d <- decide_correction(0.2, c(phase1 = 0.1))
test("ordinary improvement is corrected", stopifnot(d$decision == "corrected", d$choice == "phase1"))
d <- decide_correction(0, c(phase1 = 0))
test("zero baseline keeps the uncorrected table", stopifnot(d$decision == "uncorrected_no_batch_effect", is.na(d$choice)))
d <- decide_correction(NA_real_, c(phase1 = 0.1))
test("NA baseline is not evaluable", stopifnot(d$decision == "uncorrected_not_evaluable"))
d <- decide_correction(0.2, c(phase1 = NA_real_, phase3 = 0.08))
test("NA candidate is skipped, not fatal", stopifnot(d$decision == "corrected", d$choice == "phase3", grepl("not evaluable: phase1", d$reason)))
d <- decide_correction(0.2, c(phase1 = NA_real_, phase3 = NaN))
test("no evaluable candidate keeps the uncorrected table", stopifnot(d$decision == "uncorrected_no_improvement"))
d <- decide_correction(0.2, c(phase1 = 0.3, phase3 = 0.4))
test("all candidates worse keeps the uncorrected table", stopifnot(d$decision == "uncorrected_no_improvement"))
d <- decide_correction(1e-12, c(phase1 = 5e-13))
test("tiny values record the absolute change", stopifnot(d$decision == "corrected", grepl("absolute change -5e-13", d$reason, fixed = TRUE)))
d <- decide_correction(0.2, c(phase1 = 0.1), baseline_p = 0.3)
test("baseline batch effect not significant keeps the uncorrected table",
     stopifnot(d$decision == "uncorrected_no_batch_effect", grepl("not significant", d$reason)))
d <- decide_correction(0.2, c(phase1 = 0.1), baseline_p = 0.01)
test("significant baseline batch effect proceeds", stopifnot(d$decision == "corrected"))
d <- decide_correction(0.2, c(phase1 = 0.1, phase3 = 0.1))
test("tie keeps the earlier candidate", stopifnot(d$choice == "phase1"))
d <- decide_correction(0.2, c(phase1 = Inf, phase3 = 0.15))
test("infinite candidate never selected", stopifnot(d$choice == "phase3"))

# ---------------------------------------------------------------- run_auto_phase with injected fits
# The real control flow runs; fitting and PERMANOVA return scripted values in call order.
run_case <- function(r2, covariate_r2 = c(0.1, 0.09, 0.09), phase3_fails = FALSE) {
  calls <- 0L
  evaluate_batch_effect <<- function(data, meta, batch_col, covariates) {
    calls <<- calls + 1L
    list(batch_r2 = r2[[calls]], batch_p = 0.01, covariate_r2 = covariate_r2[[calls]], full_result = NULL)
  }
  run_phase <<- function(phase_num, tab, ...) {
    if (phase_num == 3 && phase3_fails) return(list(success = FALSE, error = "injected"))
    list(success = TRUE, data = tab + phase_num, method_matrix = matrix(phase_num), elapsed = 1)
  }
  select_optimal_reference_batch <<- function(...) "b1"
  assign("Sys.sleep", function(...) invisible(NULL), envir = globalenv())  # masks base only for this test
  tab <- matrix(1L, 20, 3, dimnames = list(sprintf("s%02d", 1:20), c("a", "b", "c")))
  meta <- data.frame(batch = rep(c("b1", "b2"), 10), age = 1:20, row.names = rownames(tab))
  capture.output(res <- run_auto_phase(tab, meta, "batch", "age", tempdir(), 0.25))
  res
}
r <- run_case(c(0.2, 0.1))
test("situation 1 ordinary stop at phase 1", stopifnot(r$decision == "corrected", r$phase_used == 1))
r <- run_case(c(0, 0, 0))
test("situation 2 zero baseline: no fit, uncorrected", stopifnot(r$decision == "uncorrected_no_batch_effect", r$phase_used == 0, all(r$data == 1L)))
r <- run_case(c(NA, 0.1, 0.08))
test("situation 3 NA baseline: no stop, uncorrected", stopifnot(r$decision == "uncorrected_not_evaluable"))
r <- run_case(c(0.2, NA, 0.08))
test("situation 4 NA phase 1: phase 3 chosen", stopifnot(r$decision == "corrected", r$phase_used == 3))
r <- run_case(c(0.2, 0.1), covariate_r2 = c(0.1, NA))
test("situation 5 NA covariate R2: no stop", stopifnot(r$decision == "corrected", r$phase_used == 1))
r <- run_case(c(1e-12, 5e-13))
test("situation 6 tiny values: corrected with absolute change", stopifnot(r$decision == "corrected", grepl("absolute change", r$decision_reason)))
r <- run_case(c(0.2, 0.3, 0.4))
test("situation 8 both worse: uncorrected, not 'success with the less bad'", stopifnot(r$decision == "uncorrected_no_improvement", r$phase_used == 0))
r <- run_case(c(0.2, 0.3), phase3_fails = TRUE)
test("phase 3 failed and phase 1 worse: uncorrected", stopifnot(r$decision == "uncorrected_no_improvement"))
r <- run_case(c(0.2, 0.18, 0.12))
test("phase 3 lower than phase 1 is chosen", stopifnot(r$decision == "corrected", r$phase_used == 3))

# ---------------------------------------------------------------- one PERMANOVA mode for every candidate
# Baseline in margin mode, phase 1 fell back to terms: compared as is, phase 1 (0.05) would stop the run.
# Re-evaluated under one mode, phase 1 is worse (0.25), phase 3 runs and wins (0.10).
mode_calls <- 0L
evaluate_batch_effect <<- function(data, meta, batch_col, covariates, mode = NULL) {
  if (identical(mode, "terms_batch_last")) {
    v <- if (all(data == 1L)) 0.2 else if (all(data == 2L)) 0.25 else 0.1
    return(list(batch_r2 = v, batch_p = 0.01, covariate_r2 = 0.1, mode = "terms_batch_last", full_result = NULL))
  }
  mode_calls <<- mode_calls + 1L
  list(batch_r2 = c(0.2, 0.05, 0.3)[mode_calls], batch_p = 0.01, covariate_r2 = 0.1,
       mode = c("margin", "terms_batch_last", "margin")[mode_calls], full_result = NULL)
}
tab40 <- matrix(1L, 20, 3, dimnames = list(sprintf("s%02d", 1:20), c("a", "b", "c")))
meta40 <- data.frame(batch = rep(c("b1", "b2"), 10), age = 1:20, row.names = rownames(tab40))
capture.output(r <- run_auto_phase(tab40, meta40, "batch", "age", tempdir(), 0.25))
test("mixed modes: all candidates re-evaluated under one mode before choosing",
     stopifnot(r$decision == "corrected", r$phase_used == 3, r$final_r2 == 0.1,
               grepl("re-evaluated with terms_batch_last", r$decision_reason)))

# ---------------------------------------------------------------- correction cohort: controls and mixed types
set.seed(1)
types <- rep(c("Primary Tumor", "Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"), 8)
batch <- rep(c("b1", "b2"), each = 16)
ids <- sprintf("S%02d", seq_along(types))
counts <- matrix(rpois(length(ids) * 3, 20) + 1L, nrow = 3, dimnames = list(NULL, ids))
otu <- data.frame(clade_name = c("k__Bacteria|g__A|s__A_one", "k__Bacteria|g__B|s__B_one", "k__Bacteria|g__C|s__C_one"),
                  counts, check.names = FALSE)
meta <- data.frame(patient = ids, type = types, batch = batch, age = rep(c(50, 61, 72, 55), 8))
o <- file.path(out, "otu.csv"); m <- file.path(out, "meta.tsv")
write.csv(otu, o, row.names = FALSE); write_tsv(meta, m)
cohort <- function(name, ..., meta_path = m) {
  a <- modifyList(list(taxon_column = "clade_name", sample_id_column = "patient", batch_column = "batch",
                       type_column = "type", tumor_only = FALSE, tumor_value = "Tumor",
                       control_values = NULL, include_controls = FALSE,
                       prefix = file.path(out, name, "fixture")), list(...))
  dir.create(dirname(a$prefix), recursive = TRUE, showWarnings = FALSE)
  assign("args", a, envir = globalenv())  # load_and_prepare_data reads the script's global args
  p <- prepare_batch_cohort(load_and_prepare_data(o, meta_path), a, "age")
  p$ledger <- read.delim(file.path(dirname(a$prefix), "batch_correction_sample_status.tsv"))
  p
}
p <- cohort("controls_default", control_values = "Blood Derived Normal")
test("controls excluded by default with a reason",
     stopifnot(!any(p$meta$type == "Blood Derived Normal"), nrow(p$meta) == 24L,
               all(grepl("^decontamination control", p$ledger$reason[p$ledger$status == "excluded"])),
               sum(p$ledger$status == "excluded") == 8L))
test("mixed tumour/normal recorded, not stopped",
     stopifnot(grepl("sample types mixed (Primary Tumor, Solid Tissue Normal)", p$sample_type_note, fixed = TRUE),
               all(p$ledger$reason[p$ledger$status == "included"] == p$sample_type_note)))
test("mixed types, both batches hold both: sample type added as a protected covariate",
     stopifnot(grepl("sample type added as a protected covariate", p$sample_type_note, fixed = TRUE),
               identical(p$covariates, c("age", "sample_type_protected")), is.factor(p$meta$sample_type_protected),
               nlevels(p$meta$sample_type_protected) == 2L))
m3 <- file.path(out, "meta_tissue_cov.tsv"); meta3 <- meta; meta3$tissue <- meta3$type; write_tsv(meta3, m3)
a3 <- list(taxon_column = "clade_name", sample_id_column = "patient", batch_column = "batch", type_column = "type",
           tumor_only = FALSE, tumor_value = "Tumor", control_values = "Blood Derived Normal", include_controls = FALSE,
           prefix = file.path(out, "tissue_cov", "fixture"))
dir.create(dirname(a3$prefix), recursive = TRUE, showWarnings = FALSE); assign("args", a3, envir = globalenv())
p3 <- prepare_batch_cohort(load_and_prepare_data(o, m3), a3, c("tissue", "age"))
test("a covariate identical to sample type: no duplicate added, recorded",
     stopifnot(identical(p3$covariates, c("tissue", "age")), grepl("already protected by covariate tissue", p3$sample_type_note)))
p <- cohort("depth_default", control_values = "Blood Derived Normal")
test("depth by batch recorded; log_depth not added by default",
     stopifnot(grepl("highest/lowest median [0-9.]+x$", p$depth_summary), !"log_depth" %in% p$covariates,
               !"log_depth" %in% names(p$meta)))
p <- cohort("depth_on", control_values = "Blood Derived Normal", depth_covariate = TRUE)
test("--depth_covariate adds log10 total count as a covariate",
     stopifnot(identical(tail(p$covariates, 1), "log_depth"),
               isTRUE(all.equal(p$meta$log_depth, log10(rowSums(p$otu)), check.attributes = FALSE))))
p <- cohort("no_protect", control_values = "Blood Derived Normal", no_protect_sample_type = TRUE)
test("protection turned off: recorded, covariates unchanged",
     stopifnot(grepl("not a protected covariate (protection turned off)", p$sample_type_note, fixed = TRUE), identical(p$covariates, "age")))
m2 <- file.path(out, "meta_single.tsv"); meta2 <- meta; meta2$type[meta2$batch == "b2" & meta2$type != "Blood Derived Normal"] <- "Primary Tumor"
write_tsv(meta2, m2)
p <- cohort("single_type_batch", control_values = "Blood Derived Normal", meta_path = m2)
test("a batch with one sample type: not added, recorded with the batch, run continues",
     stopifnot(grepl("not a protected covariate (batch b2 holds one sample type)", p$sample_type_note, fixed = TRUE), identical(p$covariates, "age")))
p <- cohort("controls_kept", control_values = "Blood Derived Normal", include_controls = TRUE)
test("--include_controls keeps controls", stopifnot(nrow(p$meta) == 32L, any(p$meta$type == "Blood Derived Normal")))
p <- cohort("no_controls_declared")
test("no declared controls: nothing excluded", stopifnot(nrow(p$meta) == 32L))
p <- cohort("tcga_tumour_only", control_values = "Blood Derived Normal", tumor_only = TRUE, tumor_value = "Primary Tumor")
test("TCGA tumour label selects tumours only, no note", stopifnot(nrow(p$meta) == 16L, all(p$meta$type == "Primary Tumor"), is.na(p$sample_type_note)))
err <- tryCatch({cohort("tumour_label_mismatch", tumor_only = TRUE, tumor_value = "Tumor"); NA_character_},
                error = function(e) conditionMessage(e))
test("wrong tumour label stops with 'No eligible samples'", stopifnot(grepl("No eligible samples", err)))

# ---------------------------------------------------------------- near-complete confounding note
mc <- data.frame(batch = rep(c("b1", "b2"), each = 20), type = c(rep("T", 19), "N", "T", rep("N", 19)),
                 age = c(rnorm(20, 40, 1), rnorm(20, 70, 1)), sex = rep(c("F", "M"), 20))
cn <- confounding_note(mc, "batch", c("type", "age", "sex"))
test("near-complete confounding is recorded for type and age, not for sex",
     stopifnot(grepl("type (Cramer's V 0.90)", cn, fixed = TRUE), grepl("age (eta squared", cn, fixed = TRUE), !grepl("sex", cn)))
test("no overlap -> NA", stopifnot(is.na(confounding_note(mc, "batch", "sex"))))

# ---------------------------------------------------------------- effect sizes for the sample-type comparison
ef <- file.path(out, "effects.tsv")
raw <- matrix(c(8, 2, 8, 2, 2, 8, 2, 8), 4, 2, byrow = TRUE, dimnames = list(c("s1", "s2", "s3", "s4"), c("tA", "tB")))
mt <- data.frame(type = c("Primary Tumor", "Primary Tumor", "Solid Tissue Normal", "Solid Tissue Normal"), row.names = rownames(raw))
write_effect_sizes(list(raw = raw, ConQuR_tuned = raw), mt, "type", "Primary Tumor", ef)
e <- read.delim(ef)
test("effect sizes: tumour value is group A, log2 ratio of mean relative abundance",
     stopifnot(all(e$status == "computed"), e$group_a[1] == "Primary Tumor", abs(e$log2_ratio_before[e$taxon == "tA"] - 2) < 1e-12,
               identical(e$log2_ratio_before, e$log2_ratio_after)))
write_effect_sizes(list(raw = raw, ConQuR_tuned = raw), data.frame(type = rep("Primary Tumor", 4), row.names = rownames(raw)), "type", "Primary Tumor", ef)
test("effect sizes: one sample type gives a reason, not an error", stopifnot(grepl("^not computed: needs exactly two", read.delim(ef)$status)))

# ---------------------------------------------------------------- effect sizes for a user-named column
# Tumour-only: sample type has one value, the named column two (s4 missing, left out of this comparison only).
mn <- data.frame(type = rep("Primary Tumor", 4), stage = c("II", "II", "I", ""), age = c(50, 60, 70, 80), row.names = rownames(raw))
write_effect_sizes(list(raw = raw, ConQuR_tuned = raw), mn, "type", "Primary Tumor", ef)
write_named_effect_sizes(list(raw = raw, ConQuR_tuned = raw), mn, "stage", c("age", "stage"), "type", ef)
e <- read.delim(ef, colClasses = "character")
st <- e[e$group_column == "stage", ]
test("named column: two values computed and appended after the sample-type row, protected recorded",
     stopifnot(nrow(e) == 3L, grepl("^not computed: needs exactly two sample types", e$status[1]), is.na(e$protected_covariate[1]),
               nrow(st) == 2L, all(st$status == "computed"), all(st$protected_covariate == "yes"), st$group_a[1] == "I", st$group_b[1] == "II",
               abs(as.numeric(st$log2_ratio_before[st$taxon == "tA"]) - log2((0.2) / 0.8)) < 1e-12))
write_effect_sizes(list(raw = raw, ConQuR_tuned = raw), mn, "type", "Primary Tumor", ef)
write_named_effect_sizes(list(raw = raw, ConQuR_tuned = raw), mn, "stage", "age", "type", ef)
test("named column not among the covariates is recorded as not protected",
     stopifnot(all(read.delim(ef, colClasses = "character")$protected_covariate[-1] == "no")))
write_effect_sizes(list(raw = raw, ConQuR_tuned = raw), mn, "type", "Primary Tumor", ef)
write_named_effect_sizes(list(raw = raw, ConQuR_tuned = raw), mn, "grade", "age", "type", ef)
e <- read.delim(ef, colClasses = "character")
test("named column missing from the metadata gives a reason, not an error",
     stopifnot(nrow(e) == 2L, e$group_column[2] == "grade", e$status[2] == "not computed: column not found in the metadata"))
mn$stage <- c("II", "II", "NA", "")
write_effect_sizes(list(raw = raw, ConQuR_tuned = raw), mn, "type", "Primary Tumor", ef)
write_named_effect_sizes(list(raw = raw, ConQuR_tuned = raw), mn, "stage", "age", "type", ef)
test("named column with one non-missing value gives a reason, not an error",
     stopifnot(grepl("^not computed: needs exactly two non-missing values among the analysis samples \\(found 1\\)",
                     read.delim(ef, colClasses = "character")$status[2])))

# write_batch_records needs vegan for its version field: run it only where the batch environment exists.
if (requireNamespace("vegan", quietly = TRUE)) {
  a <- list(prefix = file.path(out, "records", "fixture"), seed = 1, threads = 1, method = "tune",
            batch_column = "batch", control_values = "Blood Derived Normal", include_controls = FALSE)
  dir.create(dirname(a$prefix), recursive = TRUE, showWarnings = FALSE)
  tab <- matrix(1, 3, 2, dimnames = list(c("s1", "s2", "s3"), c("t1", "t2")))
  write_batch_records(list(success = TRUE, decision = "uncorrected_no_batch_effect", decision_reason = "no measurable batch effect",
                           baseline_r2 = 0, final_r2 = 0, sample_type_note = "sample types mixed (A, B); sample type is not a protected covariate"),
                      tab, c("batch", "age"), a)
  q <- read.delim(file.path(dirname(a$prefix), "batch_correction_quality.tsv"))
  s <- read.delim(file.path(dirname(a$prefix), "batch_correction_run_settings.tsv"))
  mt <- read.delim(file.path(dirname(a$prefix), "batch_correction_methods.tsv"))
  test("decision outcome recorded, not labelled corrected or failed",
       stopifnot(q$execution_status == "uncorrected_no_batch_effect", mt$selection == "uncorrected_no_batch_effect",
                 grepl("mixed", q$sample_type_note), s$value[s$key == "decontamination_controls"] == "excluded: Blood Derived Normal"))
  # Within-patient permutations (needs vegan and permute)
  suppressPackageStartupMessages(library(vegan))
  set.seed(2)
  np <- 20; ids <- c(sprintf("P%02d_T", 1:np), sprintf("P%02d_N", 1:np))
  pm <- data.frame(type = factor(rep(c("T", "N"), each = np)), age = rep(round(rnorm(np, 60, 8)), 2), row.names = ids)
  cnt <- matrix(rpois(40 * 6, 30) + 1L, 40, 6, dimnames = list(ids, paste0("t", 1:6)))
  cnt[pm$type == "T", 1:2] <- cnt[pm$type == "T", 1:2] * 3L
  patient <- setNames(substr(ids, 1, 3), ids)
  pp <- patient_restricted_p(cnt, pm, reformulate(c("type", "age")), "bray", FALSE, patient)
  test("within-patient: type varies within patients -> within-patient p; age constant within patients -> NA",
       stopifnot(pp$note == "computed", pp$n_patients == 20L, pp$n_multi == 20L, is.finite(pp$p[["type"]]), pp$p[["type"]] <= 0.05,
                 is.na(pp$p[["age"]])))
  pp1 <- patient_restricted_p(cnt[1:20, ], pm[1:20, , drop = FALSE], reformulate("age"), "bray", FALSE, patient)
  test("within-patient: one sample per patient -> not computed, with a reason", stopifnot(grepl("no patient has more than one sample", pp1$note)))
  # Same centre, different spread: the dispersion check must flag it.
  set.seed(3)
  base <- c(40, 30, 20, 10, 5, 5)
  tight <- t(replicate(20, rpois(6, base)))
  wide <- t(replicate(20, rpois(6, base * exp(rnorm(6, 0, 0.9)))))
  dm <- rbind(tight, wide) + 1L; rownames(dm) <- sprintf("d%02d", 1:40); colnames(dm) <- paste0("t", 1:6)
  dmeta <- data.frame(group = rep(c("tight", "wide"), each = 20), row.names = rownames(dm))
  dc <- dispersion_check(list(raw = dm), dmeta, "group")
  test("dispersion check: same centre, different spread -> betadisper p < 0.05",
       stopifnot(dc$Status == "computed", dc$p < 0.05, dc$F > 1))
  test("within-patient: no patient declared -> NULL", stopifnot(is.null(patient_restricted_p(cnt, pm, reformulate("type"), "bray", FALSE, NULL))))
} else cat("SKIP write_batch_records (vegan not installed here; run on the batch environment)\n")
cat("SUMMARY", passed, "tests passed\n")
