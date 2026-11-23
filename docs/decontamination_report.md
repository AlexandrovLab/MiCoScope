# Decontamination evidence in the default cohort report

The default HTML receives an explicit report bundle from the same Decontamination task. It does not search a published results directory or reuse past analyses. Reporting does not change count filtering, decontam scores, recurrence decisions or downstream batch correction.

## Reading the panels

- **Count flow:** input table → declared taxon removal → selected samples → abundance filter → prevalence filter → decontam. The cohort-selection loss comes from excluded samples; it is not a contaminant call. Nonzero taxa and table rows are reported separately.
- **Units:** in the full pipeline the numbers retain Bracken-estimated count units after upstream selection/masking. They are not counts of unique sequenced reads. For a supplied standalone table the original count semantics are not independently verified.
- **Rules:** count cells below the configured minimum become zero. Prevalence uses all selected positive and control samples, including retained zero-total samples. Taxa failing this filter are removed from the output, not merely skipped in testing.
- **Batch evidence:** sample eligibility, absent taxa, unevaluable tests and evaluated tests are distinct. The report includes detection numerators and sample denominators by taxon, batch, sample type and role at the actual test input. Whole-matrix cell detection denominators change when taxon rows are removed.
- **Controls:** labels match metadata exactly. `surrogate_biological` does not mean an experimental blank. Blood/normal tissue can have biological signal, and its comparison with tumor does not prove reagent contamination.
- **Blacklist:** exact normalized species names are flags only. The file SHA-256 identifies its content version. A missing list, an empty list and a genus-rank run are separate not-applicable states. Pooled stage shares use listed-taxon counts / all-taxon counts at that stage. Original sample shares are after declared taxon removal and before abundance/prevalence, not before all filtering.
- **Rank:** full-pipeline decontamination is species-first. Downstream genus counts are sums from retained species; a standalone genus test is a different route. Existing roll-up exclusions remain in the sample/status outputs.

## Why these rules, and what they cannot show

The existing prevalence method compares detection in configured controls and positives. The configured threshold determines per-batch calls; requiring calls in multiple batches seeks recurring evidence but can miss single-batch contamination. This reporting change does not establish an optimal threshold or validate the control choice. Defaults remain unchanged. The score histogram is descriptive, not an FDR estimate; no automatic multiple-testing correction is applied.

Setting sensitivity, depth-matched repeated resampling, alternative-mode evaluation (including `isNotContaminant`) and false-discovery calibration are not performed by this reporter. They require separate study design and adequate controls. The report labels them not assessed. Retention means the specified rules did not remove a taxon, not that it is known to be uncontaminated. A valid-empty table is a valid processing result, not evidence that a sample contains no microbes.

## Evidence and integrity

Each `*.decontam_report/` contains stage_flow.tsv, taxon_ledger.tsv, batch_type_summary.tsv, taxon_detection.tsv, the existing sample/batch/taxon decision tables, final count matrix, and manifest.tsv. The manifest records software/settings/input identities, blacklist identity and each evidence file's hash. The Python renderer checks hashes, cohort denominators, stage losses, source ledger decisions and final totals before producing HTML. A supplied incomplete/inconsistent bundle fails rendering; a legacy invocation with no bundle explicitly reports evidence unavailable.

This covers what a single run reports. It does not cover repeated stability tests, alternative mode assessment, batch-correction quality reporting, full cohort validation or release/deployment.
