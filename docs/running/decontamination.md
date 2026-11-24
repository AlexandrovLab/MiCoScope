# Decontamination

Batch-aware, prevalence-based contamination detection with the `decontam` R package
(`scripts/decontamination.R`). It is reached three ways: as part of a full run
(`--run_decontam true`), from an existing consensus table
(`--start_from decontam --consensus_otu_table <path>`), or not at all (`--run_decontam false`).
**`--run_decontam` defaults to `true`**, so a full run must either supply the parameters below
and `--metadata_file`, or pass `--run_decontam false`; otherwise it stops at launch.

**Four parameters are required whenever decontamination runs, and have no default.**
Naming which samples are the real material and which stand in for blanks is a scientific
choice, not something the pipeline should guess:

| Parameter | Meaning |
|---|---|
| `--taxonomic_rank` | The rank of the primary input table. A full run requires **`species`**, tests species and also sums retained species to genus (no second genus test). Full-run `genus` is rejected. With `--start_from decontam`, use `species` or `genus` to match the supplied table |
| `--control_role` | `experimental_blank` (true blanks) or `surrogate_biological` (a sample class standing in for blanks) |
| `--positive_values` | the sample types that are real material, i.e. **not** controls |
| `--control_values` | the sample types used as negative controls |

A run that omits any of them stops in about a second with a named error rather than
guessing. `--type_column` says which metadata column holds the role and defaults to
`sample_type`.

For TCGA, one flag supplies the whole contract — the manuscript's model, with blood-derived
normals as surrogate negative controls:

```bash
nextflow run main.nf -profile tscc --sample samples.csv \
  --metadata_file metadata.tsv --run_decontam true \
  --cohort_preset tcga --taxonomic_rank species
```

`--cohort_preset tcga` sets `--type_column "Sample Type_x"`,
`--positive_values "Primary Tumor,Solid Tissue Normal"` and
`--control_values "Blood Derived Normal"`. Passing any of them a conflicting value is an
error rather than a silent override. Note that the positives include Solid Tissue Normal:
it is real material, not a blank.

For a cohort where `Normal` means biological normal tissue (not a laboratory blank),
state that it is a surrogate control. This is a declared analysis choice, not proof that
every organism enriched in that tissue is reagent contamination:

```bash
nextflow run main.nf -profile tscc --sample samples.csv \
  --metadata_file metadata.tsv --run_decontam true --taxonomic_rank species \
  --type_column Type --control_role surrogate_biological \
  --positive_values Tumor --control_values Normal
```

For true extraction/reagent blanks instead, use their actual metadata label, for example
`--type_column Type --control_role experimental_blank --positive_values Tumor --control_values Blank`.
Do not relabel biological normal tissue or blood as `Blank` to make a run pass.

**A tumor-only final analysis still needs controls for prevalence decontamination.**
Keep the declared controls in the decontamination input; select tumors later for batch
correction with `--tumor_only true --batch_corr_tumor_values "Primary Tumor"`. The declared
controls themselves leave batch correction by default (`--downstream_include_controls false`).
That option affects batch correction, not decontamination or its output table. The TCGA
preset also retains `Solid Tissue Normal` as positive material at the decontamination
stage. To include only primary tumors and blood controls at that stage, omit the preset,
use explicit labels and exclude unselected types:

```bash
nextflow run main.nf -profile tscc --start_from decontam \
  --consensus_otu_table species_counts.tsv --decontam_library_sizes library_sizes.tsv \
  --metadata_file metadata.tsv --sample_id_column sampleid --taxonomic_rank species \
  --type_column "Sample Type_x" --batch_column shipment_batch \
  --control_role surrogate_biological --positive_values "Primary Tumor" \
  --control_values "Blood Derived Normal" --unselected_type_policy exclude
```

This example excludes `Solid Tissue Normal`; blood remains a control and remains in the
cleaned table. If the input has no appropriate controls, this prevalence analysis cannot
be interpreted as a tested tumor-only decontamination result. Do not invent controls or
change their scientific role to bypass batch eligibility checks.

A sample whose type matches neither list stops the run, because
`--unselected_type_policy` defaults to `error`; set it to `exclude` to drop those samples
instead.

The other columns and table layout:

| Parameter | Default | Meaning |
|---|---|---|
| `--metadata_file` | none (required) | per-sample metadata, tab- or comma-separated |
| `--sample_id_column` | `sampleid` | metadata column holding the sample IDs |
| `--batch_column` | `shipment_batch` | metadata column holding the batch; also used by batch correction |
| `--taxon_column` | `clade_name` | taxonomy column in the OTU table; every other column is a sample |

The ID, type and batch columns are checked against the metadata header at launch.

**How a contaminant is called.** Taxa are pre-filtered, then tested batch by batch with
decontam's prevalence method, then called across batches:

| Parameter | Default | Meaning |
|---|---|---|
| `--decontam_min_abundance` | `5` | a taxon's count in a sample below this is set to zero; count 5 is retained |
| `--decontam_min_prevalence` | `0.05` | after the count filter, a taxon must be present in at least 5% of selected samples; taxa below the threshold are excluded from both testing and the final table |
| `--decontam_threshold` | `0.1` | decontam's prevalence-test threshold |
| `--decontam_min_batches` | `2` | a taxon flagged in at least this many batches is a contaminant |
| `--min_total_per_batch`, `--min_positive_per_batch`, `--min_control_per_batch` | `5`, `1`, `1` | a batch with fewer samples, positives or controls than this is not tested |

`--decontam_min_prevalence` 0.05 matches the published method ("≥5 reads, detected in ≥5% of
samples" for CRC). It was 0.02 until 2026-09-22; results produced before then are not
comparable.

The prevalence denominator is the samples retained after role selection and library/zero-total
checks, including retained zero-count samples. At 40 retained samples, presence in 1 sample
is 2.5% (excluded), 2 is exactly 5% (eligible), and 3 is 7.5% (eligible). Passing this filter
does not guarantee that a taxon is testable in each batch or that it will be retained after
decontamination. See `filter_report.tsv` and `contaminants_final.tsv` for these separate decisions.
Count units come from the supplied table; Bracken estimates are not unique observed read counts.

Contaminants are removed **globally**, not zeroed per batch.

**Why these settings, and what they do not show.** The checks below were run
on TCGA CRC tumour and blood samples from two sequencing centres (821 samples; Broad 383 tumours and
357 blood, Baylor 32 and 49) with the default settings.

- *Prevalence mode.* The frequency mode needs a DNA-concentration measure for every library. Having
  concentrations does not by itself show that the frequency mode is better, so the prevalence mode is
  kept. `isNotContaminant` is not used: it is designed around negative (blank) controls, which TCGA
  does not have; it should be reconsidered for a cohort with extraction or library blanks.
- *Count ≥ 5 and detected in ≥ 5% of samples* follow the published CRC method; *p < 0.1* is
  decontam's default threshold for the prevalence test.
- *Called in at least 2 batches.* decontam's score is a one-sided p-value per taxon and batch, and
  each batch tests several hundred taxa, so some calls in one batch are expected by chance (at most
  about 37 of about 370 per batch at p < 0.1). Requiring a call in two batches guards against that:
  by chance about 4 taxa would be called in both, and 52 were. The calls are not FDR-corrected per
  batch; Baylor alone (81 samples, smallest p 0.004) has no call that passes Benjamini–Hochberg, but
  when each taxon's two batch p-values are combined (Fisher) and then corrected, all 52 pass at
  q < 0.05.
- *What the rule misses.* A contaminant strong in one batch only is kept: 62 taxa pass
  Benjamini–Hochberg in Broad but have p ≥ 0.1 in Baylor (mostly Bradyrhizobium, Cupriavidus,
  Burkholderia, Sphingomonas and Ralstonia). Batch-local masking, which would address them, is a
  separate decision and is not implemented.
- *Read depth.* Tumours carry 1.8–2.7× more reads than blood. Thinning counts to the blood depth
  within each batch and re-calling 20 times kept 46 of the 52 calls every time and newly called
  22 taxa from the same genera, so the depth gap seems to hide some contaminants rather than create
  them.
- *What a call means.* With blood as a surrogate control, a call says that a taxon is relatively
  common in blood compared with tumour in the same batch. It does not establish a reagent or kit
  source; the genera above are often reported as reagent contaminants, but that was not tested here.
  The cleaned table is not certified contamination-free.


**No retained taxa.** Valid input can leave no taxa after the declared filters or after
contaminant removal. The run writes a header-only count table retaining the selected sample
columns; `sample_validation.tsv`, the manifest and the HTML report distinguish `valid_empty`
with reason `no_taxa_after_filters` or `no_taxa_after_decontamination`. If no taxa reach the
statistical test, the manifest records `not_run_no_taxa_after_filters`. This is not evidence
that a sample is biologically free of microbes. Malformed counts, inconsistent declared
rank, missing metadata and the existing library/cohort guards still cause errors.

Species lacking a genus parent remain in the species table and in the genus-exclusion
ledger. If none of the surviving species can be summed to genus, the genus table is
header-only and the per-sample `genus_rollup_status.tsv` records
`not_applicable: no_genus_parent`; the HTML preserves that reason. An empty species table
instead produces `valid_empty: no_surviving_species` for the rollup.

**All-zero samples.** A sample whose counts are all zero in the table being decontaminated
(for example, a normal that carries none of the consensus genera) is evidence that those taxa
are absent, not a broken sample:

| Parameter | Default | Meaning |
|---|---|---|
| `--decontam_zero_total_policy` | `keep` | `keep`: an all-zero sample stays in the prevalence test, as long as its library was not empty. `drop`: every all-zero sample is excluded (the older behaviour, which could remove a batch's controls) |
| `--decontam_min_library` | `1000` | under `keep`, an all-zero sample whose total Bracken library is below this is excluded as `insufficient_library` |
| `--decontam_remove_taxa` | `Homo sapiens;Sinsheimervirus phiX174;synthetic construct` | species removed from the table **before** decontam (`;`-separated; at genus rank their genera, `Homo` and `Sinsheimervirus`). Human, the PhiX spike-in and *synthetic construct* (vector/adapter sequence, taxid 32630) are not microbes of the sample, so they are evidence neither for nor against contamination. Since 16 Nov 2026 human and synthetic construct never reach the analysis tables at all (review #21) and PhiX falls with the virus rule, so this list is a safety net for `--start_from decontam` with an external table. Per-sample reads removed: `reads_removed_pre_decontam` in `sample_validation.tsv`. `''` keeps them |
| `--known_contaminants` | `ref/known_contaminants.tsv` | June-Young's 95 known-contaminant species, a **flag only** (species rank): `known_contaminant` in `filter_report.tsv` / `contaminants_final.tsv`, and `known_contaminant_reads` / `pct_known_contaminant` per sample in `sample_validation.tsv`. Nothing is removed because of it. `null` turns it off |
| `--decontam_library_sizes` | none | the library-size table for `--start_from decontam`: pass the `bracken_library_sizes.tsv` that consensus wrote. In a full run it is supplied automatically |
| `--decontam_allow_unknown_library` | `false` | an all-zero sample whose library size is unknown stops the run unless this is `true`, in which case it is kept |

**Taxa that cannot be tested.** A taxon present in a batch but in fewer than 2 of its samples
cannot be tested there; decontam returns no answer.

| Parameter | Default | Meaning |
|---|---|---|
| `--decontam_unevaluable_policy` | **`skip`** | `skip`: record the taxon as `unevaluable` in that batch and not flagged there; the `--decontam_min_batches` rule then decides from the batches that could test it. `error`: stop the run |

Under `skip`, a taxon testable in fewer than `--decontam_min_batches` batches can never be called
a contaminant, however contaminant-like it looks. `contaminants_final.tsv` shows this per taxon
in `n_evaluable_batches` and `callable`; read "not a contaminant" as "could not be called" where
`callable` is `FALSE`.

**Outputs**, under `04_DECONTAMINATION/<run>/`, all prefixed with the run name:
`decontam_pkg_decontaminated.csv` (the cleaned table), `contaminants_by_batch.tsv`,
`contaminants_final.tsv`, a per-taxon `filter_report.tsv` carrying a `terminal_reason`,
`sample_validation.tsv` (every sample's role and whether it was used or excluded, and why),
`validated_metadata.tsv`, `batch_eligibility.tsv`, `run_manifest.tsv` (every setting and input
checksum), `decontamination_summary.txt` and `decontamination_plots/`.
A species-first run also writes `decontam_pkg_decontaminated.genus.csv` (surviving species summed
to genus) and `genus_rollup_excluded.tsv`: species with no genus directly above them (e.g.
endosymbionts, `synthetic_construct`) stay in the species table but are listed there, with their
counts, instead of being summed into a genus.

To decontaminate an existing consensus table without rerunning anything upstream:

```bash
nextflow run main.nf -profile tscc --start_from decontam --run_decontam true \
  --consensus_otu_table CONSENSUS_TAXA/bracken.metaphlan.common.genus.mpa.report.txt \
  --decontam_library_sizes CONSENSUS_TAXA/bracken_library_sizes.tsv \
  --metadata_file metadata.tsv --taxonomic_rank genus --cohort_preset tcga
```

---
[Back to the README](../../README.md)
