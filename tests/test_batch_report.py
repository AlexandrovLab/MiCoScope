#!/usr/bin/env python3
"""Batch-correction panel of the cohort report.

Synthetic bundles only. Checks: each run state is shown (corrected, not applied by decision, failed and
passed through); every number shown equals the bundle; the provisional review rule (covariate R2 down
>= 20%, or p from < 0.05 to >= 0.05) and its absence; evaluation-mode mismatch is not compared; a hash
or count mismatch is refused; text is escaped; no bundle means no panel.
"""
import csv
import hashlib
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("batch_report", ROOT / "scripts" / "batch_report.py")
BR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BR)


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerows(rows)


def bundle(root, state="corrected", after="ConQuR_tuned", age=(0.10, 0.09, 0.01, 0.02), sex=(0.05, 0.03, 0.04, 0.20),
           batch=(0.30, 0.05), modes=("margin", "margin"), note="NA", n=6, reason="NA", effects=None, selection=None,
           named=None):
    d = Path(root) / "run_species.batch_report"
    d.mkdir(parents=True, exist_ok=True)
    execution = "completed" if state == "corrected" else state
    write(d / "batch_correction_quality.tsv", [
        ["baseline_R2", "final_R2", "execution_status", "absolute_change", "correction_decision_reason", "sample_type_note", "n_samples", "n_taxa", "method"],
        [batch[0], batch[1], execution, batch[1] - batch[0], reason, note, n, 12, "tune"]])
    write(d / "batch_correction_run_settings.tsv", [["key", "value"], ["covariates", "age,sex"], ["batch_column", "batch"],
          ["decontamination_controls", "excluded: Blood Derived Normal"], ["ConQuR_version", "2.0"], ["vegan_version", "2.6.10"]]
          + ([["selection_distance", selection]] if selection else []))
    write(d / "batch_correction_methods.tsv", [["prevalence_bin", "selection"], ["all", "original_selected"]])
    write(d / "batch_correction_sample_status.tsv", [["sample", "status", "reason"]] + [[f"S{i}", "included", "eligible"] for i in range(n)]
          + [["B1", "excluded", "decontamination control (pass --include_controls to analyse)"]])
    perm = [["Term", "Df", "R2", "Pr(>F)", "Evaluation", "Distance", "N", "Model", "Source", "Transform"]]
    for src, mode, k in (("raw", modes[0], 0), (after, modes[1], 1)):
        if src is None:
            continue
        perm += [["batch", 1, batch[k], 0.001, mode, "bray", n, "~batch + age + sex", src, "counts"],
                 ["age", 1, age[k], age[2 + k], mode, "bray", n, "~batch + age + sex", src, "counts"],
                 ["sex", 1, sex[k], sex[2 + k], mode, "bray", n, "~batch + age + sex", src, "counts"],
                 ["batch", 1, batch[k] / 2, 0.01, mode, "euclidean", n, "~batch + age + sex", src, "RCLR"]]
    write(d / "permanova_summary.tsv", perm)
    eff_head = ["taxon", "group_column", "group_a", "group_b", "log2_ratio_before", "log2_ratio_after", "status"]
    if effects is None:
        write(d / "batch_correction_effect_sizes.tsv", [eff_head, ["NA", "type", "NA", "NA", "NA", "NA", "not computed: needs exactly two sample types (found 1)"]])
    else:
        write(d / "batch_correction_effect_sizes.tsv", [eff_head] + [[t, "type", "Tumor", "Normal", x, y, "computed"] for t, x, y in effects])
    if named is not None:       # (column, protected yes/no, [(taxon, before, after)] or a "not computed: ..." status)
        col, prot, rows = named
        old = list(csv.reader((d / "batch_correction_effect_sizes.tsv").open(), delimiter="\t"))
        new = [old[0] + ["protected_covariate"]] + [r + ["NA"] for r in old[1:]]
        new += ([["NA", col, "NA", "NA", "NA", "NA", rows, prot]] if isinstance(rows, str) else
                [[t, col, "II", "I", x, y, "computed", prot] for t, x, y in rows])
        write(d / "batch_correction_effect_sizes.tsv", new)
    write(d / "batch_correction_dispersion.tsv", [["Source", "Group", "F", "p", "N", "Distance", "Status"],
          ["raw", "batch", 12.5, 0.001, n, "bray", "computed"], [after, "batch", 0.8, 0.41, n, "bray", "computed"]])
    write(d / "BATCH_CORRECTION_STATUS.tsv", [["status", "reason", "n_samples"], [state, reason, n]])
    rows = [["key", "value"], ["report_schema", "1"], ["rank", "species"]]
    rows += [[f"report_sha256.{f}", hashlib.sha256((d / f).read_bytes()).hexdigest()] for f in BR.FILES]
    write(d / "manifest.tsv", rows)
    return d


class BatchPanel(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_corrected_numbers_and_review_rule(self):
        html = BR.render([bundle(self.tmp)])
        self.assertIn("Corrected with ConQuR (tune)", html)
        self.assertIn("<td>0.3000</td><td>0.0500</td><td>-0.2500</td>", html)       # batch before/after/change
        rows, after = BR.term_rows(BR.read(Path(self.tmp) / "run_species.batch_report" / "permanova_summary.tsv"), "batch", ["age", "sex"])
        flags = {r["term"]: r["review"] for r in rows}
        self.assertEqual(after, "ConQuR_tuned")
        # Batch R2 0.30 -> 0.05, so covariates are judged on their share of the non-batch variation:
        self.assertTrue(flags["age"])           # R2 0.10 -> 0.09, share 0.143 -> 0.095 (a 34% fall)
        self.assertTrue(flags["sex"])           # share 0.071 -> 0.032 and p 0.04 -> 0.20
        self.assertIn("Review needed: age, sex.", html)
        self.assertIn("0.1429 → 0.0947", html)  # the share column
        self.assertIn("excluded: Blood Derived Normal", html)
        self.assertIn("RCLR", html)              # transform table folded in a details block
        self.assertIn("Not assessed in this run", html)

    def test_share_of_non_batch_variation(self):
        # A rise that only reflects the removed batch variation is not flagged; a real rise in share is.
        d = bundle(self.tmp, batch=(0.72, 0.0), age=(0.041, 0.156, 0.30, 0.01), sex=(0.010, 0.060, 0.30, 0.01))
        rows, _ = BR.term_rows(BR.read(d / "permanova_summary.tsv"), "batch", ["age", "sex"])
        flags = {r["term"]: (r["review"], r["note"]) for r in rows}
        self.assertFalse(flags["age"][0])   # share 0.146 -> 0.156: the rise is the batch removal
        self.assertTrue(flags["sex"][0]); self.assertIn("share of non-batch variation rose", flags["sex"][1])
        d2 = bundle(tempfile.mkdtemp(), batch=(0.72, 0.0), age=(-0.002, 0.03, 0.9, 0.5), sex=(0.05, 0.01, 0.6, 0.7))
        rows2, _ = BR.term_rows(BR.read(d2 / "permanova_summary.tsv"), "batch", ["age", "sex"])
        self.assertFalse(any(r["review"] for r in rows2 if r["role"] == "covariate"))   # not significant either side

    def test_significance_loss_alone_is_flagged(self):
        d = bundle(self.tmp, age=(0.10, 0.095, 0.04, 0.06), sex=(0.05, 0.049, 0.30, 0.30), batch=(0.30, 0.30))
        rows, _ = BR.term_rows(BR.read(d / "permanova_summary.tsv"), "batch", ["age", "sex"])
        flags = {r["term"]: r["review"] for r in rows}
        self.assertTrue(flags["age"])            # 5% fall, but p 0.04 -> 0.06
        self.assertFalse(flags["sex"])

    def test_uncorrected_by_decision(self):
        d = bundle(self.tmp, state="uncorrected_no_batch_effect", after="UNCORRECTED_by_decision", batch=(0.0, 0.0),
                   age=(0.10, 0.10, 0.01, 0.01), sex=(0.05, 0.05, 0.04, 0.04), reason="no measurable batch effect before correction")
        html = BR.render([d])
        self.assertIn("Correction not applied (uncorrected_no_batch_effect)", html)
        self.assertIn("The input table is used downstream", html)
        self.assertIn("No covariate met the review rule.", html)

    def test_correction_decision_basis_is_shown(self):
        # The panel states the relative-abundance decision rule and no longer lists
        # relative-abundance Bray as not assessed; the R2 table itself stays on counts.
        rule = "bray on relative abundance; correction starts only when the baseline batch p < 0.05"
        html = BR.render([bundle(self.tmp, selection=rule)])
        self.assertIn(f"Correction decision: {rule}.", html.replace("&lt;", "<"))
        self.assertNotIn("Relative-abundance Bray distances", html)
        self.assertIn("Compositional metrics (for example Aitchison distance) and stratification by library size.", html)
        self.assertNotIn("Correction decision:", BR.render([bundle(tempfile.mkdtemp())]))   # older bundles: no line

    def test_decision_metric_is_shown(self):
        # G2 (2026-11-20): the relative-abundance numbers the correction decision used are shown.
        d = bundle(self.tmp)
        q = BR.read(d / "batch_correction_quality.tsv")[0]
        cols = list(q.keys()) + ["baseline_p"]
        write(d / "batch_correction_quality.tsv", [cols, [q[k] for k in q] + ["0.001"]])
        m = [["key", "value"], ["report_schema", "1"], ["rank", "species"]] + [[f"report_sha256.{f}", hashlib.sha256((d / f).read_bytes()).hexdigest()] for f in BR.FILES]
        write(d / "manifest.tsv", m)
        html = BR.render([d])
        self.assertIn("Decision metric (Bray–Curtis on relative abundances): batch R2 0.3000 (p 0.001) before; 0.0500 after correction.", html)
        self.assertNotIn("Decision metric", BR.render([bundle(tempfile.mkdtemp())]))   # older bundles: no line

    def test_depth_by_batch_is_shown(self):
        # Depth by batch is shown with guidance; when depth was added, the panel says so.
        d = bundle(self.tmp)
        rows = BR.read(d / "batch_correction_run_settings.tsv")
        def rewrite(extra):
            write(d / "batch_correction_run_settings.tsv", [["key", "value"]] + [[r["key"], r["value"]] for r in rows] + extra)
            m = [["key", "value"], ["report_schema", "1"], ["rank", "species"]] + [[f"report_sha256.{f}", hashlib.sha256((d / f).read_bytes()).hexdigest()] for f in BR.FILES]
            write(d / "manifest.tsv", m)
            return BR.render([d])
        html = rewrite([["depth_by_batch", "B1 2,000; B2 1,000; highest/lowest median 2.00x"], ["depth_covariate", "not added"]])
        self.assertIn("Sequencing depth (median total count by batch): B1 2,000; B2 1,000; highest/lowest median 2.00x.", html)
        self.assertIn("--batch_corr_depth_covariate true", html)
        html = rewrite([["depth_by_batch", "B1 2,000; B2 1,000; highest/lowest median 2.00x"], ["depth_covariate", "added (log_depth = log10 total count)"]])
        self.assertIn("Depth (log10 total count) was added as a covariate.", html)
        self.assertNotIn("--batch_corr_depth_covariate true", html)
        self.assertNotIn("Sequencing depth (median", BR.render([bundle(tempfile.mkdtemp())]))   # older bundles: no line

    def test_mixed_type_note_is_shown(self):
        note = "sample types mixed (Primary Tumor, Solid Tissue Normal); sample type is not a protected covariate"
        html = BR.render([bundle(self.tmp, note=note)])
        self.assertIn(note.replace("(", "(").replace(")", ")"), html)
        self.assertIn("the run continued as configured", html)

    def test_confounding_note_shown(self):
        d = bundle(self.tmp)
        q = BR.read(d / "batch_correction_quality.tsv")[0]
        cols = list(q.keys()) + ["confounding_note"]
        write(d / "batch_correction_quality.tsv", [cols, [q[k] for k in q] + ["batch nearly overlaps with type (Cramer's V 0.90); batch and these effects may not be separable"]])
        m = [["key", "value"], ["report_schema", "1"], ["rank", "species"]] + [[f"report_sha256.{f}", hashlib.sha256((d / f).read_bytes()).hexdigest()] for f in BR.FILES]
        write(d / "manifest.tsv", m)
        html = BR.render([d])
        self.assertIn("Review needed: batch nearly overlaps with type (Cramer&#x27;s V 0.90)", html)

    def test_protected_note_has_no_warning(self):
        note = "sample types mixed (A, B); sample type added as a protected covariate (sample_type_protected)"
        html = BR.render([bundle(self.tmp, note=note)])
        self.assertIn(note, html)
        self.assertNotIn("shrank a known tissue difference", html)

    def test_patient_permutation_shown(self):
        d = bundle(self.tmp)
        rows = BR.read(d / "permanova_summary.tsv")
        head = list(rows[0].keys()) + ["p_within_patient", "Patients", "Patients_multi", "Patient_note"]
        out = [head]
        for r in rows:
            pp = "0.004" if r["Term"] == "sex" else "NA"
            out.append([r[k] for k in rows[0].keys()] + [pp, 30, 12, "computed"])
        write(d / "permanova_summary.tsv", out)
        m = [["key", "value"], ["report_schema", "1"], ["rank", "species"]] + [[f"report_sha256.{f}", hashlib.sha256((d / f).read_bytes()).hexdigest()] for f in BR.FILES]
        write(d / "manifest.tsv", m)
        html = BR.render([d])
        self.assertIn("p before (within patient)", html)
        self.assertIn("Within-patient permutation: 30 patients, 12 with more than one sample.", html)
        self.assertNotIn("Patient-restricted permutations.", html)

    def test_no_patient_no_column(self):
        html = BR.render([bundle(self.tmp)])
        self.assertNotIn("within patient", html)
        self.assertIn("Patient-restricted permutations.", html)

    def test_dispersion_shown(self):
        html = BR.render([bundle(self.tmp)])
        self.assertIn("Dispersion check (betadisper, Bray): batch: F 12.50 → 0.80, p 0.001 → 0.410", html)
        self.assertNotIn("Within-group dispersion", html)

    def test_mode_mismatch_not_compared(self):
        d = bundle(self.tmp, modes=("margin", "terms_batch_last"))
        rows, _ = BR.term_rows(BR.read(d / "permanova_summary.tsv"), "batch", ["age", "sex"])
        self.assertTrue(all(r["change"] is None and "not directly comparable" in r["note"] for r in rows))
        self.assertFalse(any(r["review"] for r in rows))

    def test_not_computable_is_noted_not_flagged(self):
        d = bundle(self.tmp, sex=(0.05, "NA", 0.04, "NA"))
        rows, _ = BR.term_rows(BR.read(d / "permanova_summary.tsv"), "batch", ["age", "sex"])
        sex = next(r for r in rows if r["term"] == "sex")
        self.assertFalse(sex["review"])          # only the two stated triggers set the flag
        self.assertIn("R2 not computable", sex["note"])

    def test_failed_passthrough_shown(self):
        d = bundle(self.tmp, state="failed", after="UNCORRECTED_passthrough", reason="ConQuR fit error")
        html = BR.render([d])
        self.assertIn("Status failed: ConQuR fit error.", html)
        self.assertIn("After-correction table: UNCORRECTED_passthrough.", html)

    def test_effect_size_shrinkage_flagged(self):
        eff = [("t1", 2.0, 1.2), ("t2", -2.0, -1.4), ("t3", 0.5, 0.1)]          # strong: t1, t2; kept 2.6/4 = 65%
        html = BR.render([bundle(self.tmp, effects=eff)])
        self.assertIn("Sample-type comparison: Tumor vs Normal (type)", html)
        self.assertIn("Review needed: 2 taxa differed at least two-fold", html)
        self.assertIn("is 65% of the value before, and 0 changed direction", html)
        self.assertNotIn("none is specified for this run", html)          # the named comparison is now assessed

    def test_effect_size_direction_change_flagged(self):
        e, _ = BR.effect_summary(BR.read(bundle(self.tmp, effects=[("t1", 1.5, -1.6), ("t2", 2.0, 2.1)]) / "batch_correction_effect_sizes.tsv"))
        self.assertEqual(e["flips"], 1)
        self.assertTrue(e["review"])

    def test_effect_size_kept_not_flagged(self):
        e, _ = BR.effect_summary(BR.read(bundle(self.tmp, effects=[("t1", 2.0, 1.9), ("t2", -1.0, -0.95)]) / "batch_correction_effect_sizes.tsv"))
        self.assertFalse(e["review"])

    def test_effect_size_not_computed_reason_shown(self):
        html = BR.render([bundle(self.tmp)])
        self.assertIn("not computed: needs exactly two sample types (found 1).", html)
        self.assertIn("none is specified for this run", html)

    def test_named_comparison_rendered_next_to_sample_type(self):
        html = BR.render([bundle(self.tmp, effects=[("t1", 2.0, 1.9)], named=("stage", "yes", [("t1", 1.5, -1.6), ("t2", 2.0, 2.1)]))])
        self.assertIn("Sample-type comparison: Tumor vs Normal (type)", html)
        self.assertIn("Named comparison: II vs I (stage)", html)
        self.assertEqual(html.count("Largest differences before correction"), 2)
        self.assertIn("Review needed: 2 taxa differed at least two-fold", html)     # one direction change in stage
        self.assertIn("No review flag: 1 taxa differed", html)                    # sample type kept
        self.assertIn("stage is a protected covariate", html)
        self.assertNotIn("not a protected covariate", html)
        self.assertIn("Clinical comparisons other than sample type and stage.", html)

    def test_named_comparison_in_tumour_only_run(self):
        html = BR.render([bundle(self.tmp, named=("stage", "no", [("t1", 2.0, 1.0)]))])
        self.assertIn("not computed: needs exactly two sample types (found 1).", html)
        self.assertIn("Named comparison: II vs I (stage)", html)
        self.assertIn("stage is not a protected covariate: correction can shrink this difference; "
                      "list it in --batch_corr_covariates to protect it.", html)
        self.assertNotIn("none is specified for this run", html)
        self.assertIn("Clinical comparisons other than stage.", html)

    def test_named_comparison_not_computed_reason_shown(self):
        html = BR.render([bundle(self.tmp, named=("stage", "no", "not computed: column not found in the metadata"))])
        self.assertIn("<h3>Named comparison (stage)</h3><p>not computed: column not found in the metadata.</p>", html)
        self.assertNotIn("protected covariate", html)      # protection matters only for a computed comparison

    def test_older_bundle_without_named_rows_unchanged(self):
        eff = [("t1", 2.0, 1.2), ("t2", -2.0, -1.4)]
        old = BR.render([bundle(self.tmp, effects=eff)])
        new = BR.render([bundle(tempfile.mkdtemp(), effects=eff, named=("stage", "yes", [("t1", 1.0, 1.0)]))])
        self.assertNotIn("protected_covariate", old)
        self.assertNotIn("Named comparison", old)
        self.assertEqual(old.split("<h3>Not assessed")[0], new.split("<h3>Named comparison")[0])

    def test_mixed_warning_shown_with_note(self):
        html = BR.render([bundle(self.tmp, note="sample types mixed (A, B); sample type is not a protected covariate")])
        self.assertIn("shrank a known tissue difference", html)

    def test_hash_mismatch_refused(self):
        d = bundle(self.tmp)
        with (d / "permanova_summary.tsv").open("a") as f:
            f.write("extra\n")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            BR.render([d])

    def test_count_mismatch_refused(self):
        d = bundle(self.tmp)
        write(d / "batch_correction_quality.tsv", [["n_samples", "n_taxa", "method", "execution_status", "sample_type_note"], [99, 12, "tune", "completed", "NA"]])
        rows = [["key", "value"], ["report_schema", "1"], ["rank", "species"]]
        rows += [[f"report_sha256.{f}", hashlib.sha256((d / f).read_bytes()).hexdigest()] for f in BR.FILES]
        write(d / "manifest.tsv", rows)
        with self.assertRaisesRegex(ValueError, "sample count disagrees"):
            BR.render([d])

    def test_text_is_escaped(self):
        html = BR.render([bundle(self.tmp, state="uncorrected_no_improvement", after="UNCORRECTED_by_decision", reason="<b>x</b>")])
        self.assertIn("&lt;b&gt;x&lt;/b&gt;", html)
        self.assertNotIn("<b>x</b>", html)

    def test_no_bundle_no_panel(self):
        self.assertEqual(BR.render([]), "")


if __name__ == "__main__":
    unittest.main()
