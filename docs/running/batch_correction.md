# Batch correction (opt-in)

`--run_batch_correction true` removes batch effects from the decontaminated table with ConQuR
(`scripts/2500703_batch_correction_normalization.r`), then normalises it and measures how much
variance batch explains before and after. It takes its input from decontamination, or from
consensus/Bracken when decontamination is off, or from an existing table with
`--start_from batch_correction --decontam_otu_table <path> --metadata_file <path>`. It uses
`--batch_column`, `--sample_id_column` and `--type_column` from the decontamination settings.

| Parameter | Default | Meaning |
|---|---|---|
| `--batch_corr_covariates` | `age_diag,sex,bmi` | comma-separated metadata columns ConQuR adjusts for. Every named column must exist in the metadata, or the run stops with a named error |
| `--batch_corr_method` | `tune` | `tune` (auto-phased Tune_ConQuR) or `vanilla` (default ConQuR settings). Nothing else is accepted: ComBat is not implemented |
| `--phase` | `auto` | only `auto` is implemented: phase 1 (16 settings), then phase 3 if phase 1 did not reduce the batch effect enough |
| `--r2_threshold` | `0.25` | phase 1 is accepted if it cuts the batch R² by at least this fraction |
| `--tumor_only` | `false` | correct tumours only. Which types count as tumours is `--batch_corr_tumor_values` |
| `--batch_corr_tumor_values` | `Tumor`; `Primary Tumor` with `--cohort_preset tcga` | comma-separated, case-insensitive. Not the same as decontamination's `--positive_values`, which means "not a control" |
| `--downstream_include_controls` | `false` | decontamination controls (`--control_values`) are left out of batch correction and the analyses after it, and listed in `batch_correction_sample_status.tsv` with the reason. `true` keeps them, for a design that deliberately analyses them |
| `--batch_corr_allow_uncorrected` | `false` | a failed correction stops the run. `true` continues with the raw counts, published as `corrected/UNCORRECTED_passthrough.tsv`, never under a corrected name |
| `--batch_corr_env` | built by Nextflow from the batch lockfile | optional pre-built prefix; see setup guide below |
| `--batch_corr_seed` | `1` | RNG seed; hold allocated CPUs and environment fixed for repeat comparisons |
| `--conqur_sha` | `ff233085…` | the ConQuR commit the script was written against; any other install is refused |

## Choose the analysis cohort explicitly

These examples start from a species table that has already been decontaminated; they do
not rerun upstream classification. Replace paths and column names with those in your files.
Here `age_diag` and `sex` must exist and satisfy the checks below; these examples do not
require a BMI column. Use a sample/library identifier for `--sample_id_column`, not a
patient identifier shared by multiple libraries.

**TCGA primary tumors only:**

```bash
nextflow run main.nf -profile tscc --start_from batch_correction \
  --run_decontam false --run_batch_correction true \
  --decontam_otu_table species_counts.cleaned.csv --metadata_file metadata.tsv \
  --taxonomic_rank species --sample_id_column sampleid --batch_column shipment_batch \
  --type_column "Sample Type_x" --tumor_only true \
  --batch_corr_tumor_values "Primary Tumor" --batch_corr_covariates age_diag,sex
```

`Tumor` does not match `Primary Tumor`; with `--cohort_preset tcga` the default is already
`Primary Tumor`. Here primary tumors are eligible for further checks;
`Solid Tissue Normal` and `Blood Derived Normal` are excluded from batch correction and
remain in the sample status ledger with the reason. Tumor selection is not applied to
any preceding decontamination step. For a separate genus table use `--taxonomic_rank genus`.

**Mixed tissues:** use the same command with `--tumor_only false`, omitting the tumor-values option.
This includes primary tumor and solid normal when present and otherwise eligible. Blood used
as a decontamination control is left out by default; add `--downstream_include_controls true`
only if the design analyses it on purpose. When several sample types remain and sample type is not
among the covariates, it is added as the protected covariate `sample_type_protected` if every batch
holds at least two types (`--batch_corr_protect_sample_type`, default `true`), unless a listed
covariate already separates the types exactly. If a batch holds only one type, or protection is
turned off, the run continues without it. Either way the choice is recorded on every included sample
in `batch_correction_sample_status.tsv`, in `batch_correction_quality.tsv` (`sample_type_note`) and in
the cohort report. Keeping tissue type in the model does not prove that biological signals were
preserved.

A missing requested column is an error. Missing values, unknown batches and all-zero count
samples can exclude samples; read `batch_correction_sample_status.tsv`. Do not add invented
BMI values or silently remove covariates to obtain a successful run. If every batch has only
one tissue type, tissue and batch effects cannot be separated; the existing within-batch
covariate check stops such a request. Near confounding and repeated patients need separate
design and evaluation work; a successful fit alone is not evidence of identifiability.

## When batches differ in sequencing depth

Whether to correct is decided on Bray–Curtis of relative abundances, and correction runs only when the
baseline batch p is below 0.05, so a depth difference alone does not trigger correction. When batches
differ both in composition and in depth, correction can also remove part of a true biological
difference. The cohort report therefore shows the median total count per batch and the ratio of the
highest to the lowest median, and leaves the choice to the user:

| Choice | Effect in a synthetic check (2026-11-19; a real batch shift plus a depth gap, 6 runs each) |
|---|---|
| Default (no depth term) | 57–76% of a known tissue difference kept; batch effect removed |
| `--batch_corr_depth_covariate true` (adds `log_depth` = log10 total count per sample) | 83–106% kept; 9–11% of the before-correction batch R² left at a 1.77× depth gap, 36–45% at 2.65× |
| Rarefying before correction (not built in) | 88–114% kept; batch effect removed; discards reads |

Depth tracks batch in such data, so protecting depth also protects part of the batch difference; the
larger the gap, the more is left. With the option on, `log_depth` also enters the near-overlap check, so a
strong depth–batch association is flagged “review needed”. Synthetic data, one design and three seeds:
treat these numbers as a guide, not a calibration.

## Over-correction: flags and sign-off, not an automatic block

Batch correction can shrink a real biological difference. The pipeline does not stop or reject a
correction on its own, because no threshold has been validated: in
a synthetic check, covariate R² caught 1 of 10 real losses, and the per-taxon effect rule caught all
10 but also flagged most runs where batch and sample type overlap. Instead, the cohort report shows
for each rank:

- batch and covariate R² and p before and after correction, with the absolute change, and each
  covariate's share of the non-batch variation (removing batch variation makes every other R² grow by
  itself, so the share is the comparable number);
- the per-taxon effect size of the sample-type comparison before and after;
- “review needed” when a covariate's share of the non-batch variation, R² ÷ (1 − batch R²), falls by
  at least 20% (for a covariate with p < 0.05 before correction) or rises by at least 20% and 0.01 (with
  p < 0.05 after), or its p moves from below 0.05 to 0.05 or above, the effect size falls by at least 20% or changes sign, or batch nearly overlaps with a
  covariate (Cramér's V or η² ≥ 0.8). These thresholds are provisional.

In a tumour-only run there is only one sample type, so that comparison is not computed. Name a metadata
column with two values, such as a clinical group, with `--batch_corr_effect_column`; its per-taxon effect
size before and after correction is shown with the same flag. If the column is not listed in
`--batch_corr_covariates`, correction does not protect it and can shrink the difference, and the report
says so.

Before results from a run are released or used in a manuscript, a named reviewer reads every
“review needed” flag and records, in the release notes or analysis log, the run, each flag, the
decision (accept the corrected table, use the uncorrected table, or change the model) and the
reviewer and date. A flag is a pointer for that review, not a pass or a fail.

Sample type is protected automatically when every batch holds at least two types; a
single-type batch cannot be protected, and the report says so. To state protection explicitly, add a
sample-type column whose name has no spaces to the metadata and list it in `--batch_corr_covariates`;
the automatic copy is then not added (the TCGA column `Sample Type_x` cannot be used in a model formula
as is).

## Comparing settings

The cohort report describes one run. The pipeline has no single command that compares filter,
decontamination, batch-correction or normalization settings. To compare settings, run each one on the
same input table and samples, with the same seed, into its own output folder, and keep the default run
as the baseline. Starting from a saved table (`--start_from decontam` or `batch_correction`) is much
cheaper than rerunning classification. Differences between runs show how sensitive the calls are to a
setting, not which result is correct.

**No taxa to correct.** Before scheduling ConQuR, the pipeline checks each rank's table.
A valid header-only table writes `empty_counts.csv` with the original sample columns,
`<rank>.batch_sample_status.tsv` and `batch_skip_status.tsv` under that rank's output folder.
The latter records correction, normalization and PERMANOVA as `skipped` with reason
`no_taxa_in_input`; the HTML displays the skip. No corrected matrix or statistical result is
invented. Nonempty tables continue through the existing correction path; malformed tables
and invalid metadata remain errors.

**What ConQuR needs**, and what happens without it:

- **At least one covariate that varies within every batch.** `--batch_corr_covariates none`
  fails loudly, by design: ConQuR builds its model from covariates plus batch and cannot run on
  batch alone. A requested covariate that is constant within any one eligible batch also stops the
  run, naming the covariate and batch: the model is never changed silently, so change
  `--batch_corr_covariates` explicitly if that covariate cannot be used.
- **At least 20 samples** for the tuned method, after these removals, each logged: samples that
  are all-zero after decontamination (listed in `batch_correction_sample_status.tsv`), samples
  with an unknown or missing batch, samples missing any covariate value, and batches with fewer
  than 5 samples.

A failure writes `corrected/BATCH_CORRECTION_STATUS.tsv` (`failed`, the reason and sample counts)
and stops the run, unless `--batch_corr_allow_uncorrected true`.

**When correction is not applied.** These are recorded outcomes, not failures, and the run continues
with the input table published as `corrected/UNCORRECTED_by_decision.tsv` (never as `ConQuR_tuned.tsv`):
`uncorrected_no_batch_effect` when the batch R² before correction is 0 or below;
`uncorrected_not_evaluable` when that R² cannot be computed; `uncorrected_no_improvement` when no
candidate with a computable batch R² lowers it. A candidate whose batch R² cannot be computed is never
selected. The status file, `batch_correction_quality.tsv` (`execution_status`,
`correction_decision_reason`, `absolute_change`) and the parameters file give the reason.

**Ranks.** Species-first runs (a full run, or `--start_from decontam` with a species table) correct
the species table and the genus table summed from it, as two independent ConQuR runs with the same
samples, covariates and seed. A supplied table (`--start_from batch_correction`, or a standalone
genus decontamination) is corrected at its own `--taxonomic_rank`.

**Outputs**, under `05_BATCH_CORRECTION/<run>_<rank>/` (e.g. `decontam_run_species/`, `decontam_run_genus/`): `corrected/ConQuR_tuned.tsv` and
`ConQuR_tuned_parameters.txt` (the file keeps this name under `--batch_corr_method vanilla`
too; its parameters file then reads `Phase used: vanilla`), `corrected/BATCH_CORRECTION_STATUS.tsv`,
`normalized/` (the raw and corrected tables under several transforms: CLR with a pseudocount, CLR
with zCompositions count-zero multiplicative replacement, and RCLR; the RCLR table followed by
Euclidean distance is a sensitivity analysis, not robust PCA (RPCA), which would add matrix completion), `pcoa_plots/`,
`permanova_summary.tsv` (batch and covariate R² for every table and transform; check this to
judge the correction), `batch_correction_sample_status.tsv`, `batch_correction_taxon_status.tsv`,
`batch_correction_quality.tsv`, `batch_correction_methods.tsv` and `batch_correction_run_settings.tsv`.

---
[Back to the README](../../README.md)

## Fresh environment and reproducible execution

See [batch setup and reproducibility](batch_setup_and_reproducibility.md) for the official bootstrap command, fixed-seed execution, quality reports and validation scope.

## ConQuR compatibility build

The batch lock pins ConQuR2.0 at upstream `ff233085cedc36a24382038067a4cece3e94e0a0`
plus the separately identified `adonis2-terms-v1` patch (conda build `r43_1`) and vegan2.6-10.
The patch replaces the four PERMANOVA scorer calls used by `Tune_ConQuR` with
`adonis2(method="bray", permutations=999, by="terms")` and named R2 extraction.
Sequential evaluation preserves the legacy selection definition; this is not a new
criterion for biological signal preservation. Other exported ConQuR functions are not patched.
The single-batch runtime guard still preserves original counts when a prevalence bin cannot be scored.
Both upstream SHA and local patch identity are checked before the batch script runs and
recorded in `batch_correction_run_settings.tsv` along with the actual vegan version.
Existing manually created environments must be rebuilt with `scripts/setup_batch_env.sh`
into a new prefix; they are not modified automatically.
