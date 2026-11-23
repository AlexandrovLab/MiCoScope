"""Behavioral contracts: valid empty != malformed; skipped correction reaches HTML."""
import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import test_decontam_zero_total as zero

ROOT = Path(__file__).resolve().parents[1]


def run_script(name, cwd, *args):
    return subprocess.run([sys.executable, str(ROOT / "scripts" / name), *map(str, args)],
                          cwd=cwd, capture_output=True, text=True)


class EmptyRouting(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        (self.d / "meta.tsv").write_text("sample\tbatch\ttype\nX\tB\tTumor\n")

    def tearDown(self):
        self.tmp.cleanup()

    def prepare(self, text, rank="genus"):
        (self.d / "input.csv").write_text(text)
        return run_script("prepare_batch_input.py", self.d, "--input", "input.csv", "--metadata", "meta.tsv",
                          "--rank", rank, "--sample-column", "sample", "--batch-column", "batch", "--type-column", "type")

    def test_empty_table_retains_sample_and_skips_without_fake_results(self):
        p = self.prepare("clade_name,X\n")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertFalse((self.d / "ready_counts.tsv").exists())
        self.assertEqual((self.d / "empty_counts.csv").read_text(), "clade_name,X\n")
        self.assertIn("X\tskipped\tno_taxa_in_input", (self.d / "genus.batch_sample_status.tsv").read_text())
        self.assertIn("permanova\tskipped", (self.d / "batch_skip_status.tsv").read_text())
        self.assertFalse((self.d / "permanova_summary.tsv").exists())
        self.check_report("skipped: no_taxa_in_input")

    def check_report(self, expected):
        p = run_script("build_qc_summary.py", self.d, "--inputs", ".", "--output", "qc.tsv")
        self.assertEqual(p.returncode, 0, p.stderr)
        p = run_script("build_cohort_report.py", self.d, "--qc-summary", "qc.tsv", "--output", "report.html")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn(expected, (self.d / "report.html").read_text())

    def test_no_genus_has_its_own_reason_and_excluded_counts(self):
        (self.d / "species.csv").write_text("clade_name,X\ns__unknown,8\n")
        p = run_script("aggregate_decontaminated_species.py", self.d, "--input", "species.csv",
                       "--output", "genus.csv", "--excluded", "excluded.tsv", "--status", "run.genus_rollup_status.tsv")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("no genus in lineage\t8\t1", (self.d / "excluded.tsv").read_text())
        (self.d / "run.sample_validation.tsv").write_text(
            "sample_id\tin_otu\tsample_type\trole\tfinal_status\nX\tTRUE\tTumor\tpositive\tincluded\n")
        self.check_report("not_applicable: no_genus_parent")

    def test_nonempty_input_passes_through_byte_for_byte(self):
        text = "clade_name,X\ng__G,8\n"
        p = self.prepare(text)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual((self.d / "ready_counts.tsv").read_text(), text)
        self.assertFalse((self.d / "genus.batch_sample_status.tsv").exists())

    def test_malformed_counts_or_taxonomy_cannot_become_a_skip(self):
        for text in ["", "clade_name,X\ng__G,-1\n", "clade_name,X\ng__G,NaN\n",
                     "clade_name,X\ns__wrong_rank,0\n", "clade_name,X\ng__G,1,2\n"]:
            with self.subTest(text=text):
                p = self.prepare(text)
                self.assertNotEqual(p.returncode, 0)
                self.assertFalse((self.d / "genus.batch_sample_status.tsv").exists())

    def test_empty_input_does_not_hide_invalid_metadata(self):
        (self.d / "meta.tsv").write_text("sample\tbatch\ttype\nY\tB\tTumor\n")
        p = self.prepare("clade_name,X\n")
        self.assertNotEqual(p.returncode, 0)


@unittest.skipUnless(zero.decontam_available(), "R/decontam/phyloseq/optparse/digest test dependencies unavailable")
class EmptyDecontamination(unittest.TestCase):
    setUp = zero.ZeroTotalBehaviour.setUp
    tearDown = zero.ZeroTotalBehaviour.tearDown
    run_r = zero.ZeroTotalBehaviour.run_r

    def test_all_filtered_writes_valid_empty_and_html_reason(self):
        p = self.run_r("--min_abundance", "1000")
        self.assertEqual(p.returncode, 0, p.stderr)
        d = Path(self.dir)
        with (d / "t.decontam_pkg_decontaminated.csv").open() as f:
            rows = list(csv.reader(f))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][1:], [f"S{i}" for i in range(1, 13)])
        with (d / "t.sample_validation.tsv").open() as f:
            statuses = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual({r["result_status"] for r in statuses}, {"valid_empty"})
        self.assertEqual({r["result_reason"] for r in statuses}, {"no_taxa_after_filters"})
        self.assertIn("not_run_no_taxa_after_filters", (d / "t.run_manifest.tsv").read_text())
        p = run_script("build_qc_summary.py", d, "--inputs", ".", "--output", "qc.tsv")
        self.assertEqual(p.returncode, 0, p.stderr)
        p = run_script("build_cohort_report.py", d, "--qc-summary", "qc.tsv", "--output", "report.html")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("no_taxa_after_filters", (d / "report.html").read_text())
        self.assertIn("completed_valid_empty_decontamination", (d / "qc.tsv").read_text())

    def test_explicit_removal_of_every_taxon_preserves_all_table_headers(self):
        p = self.run_r("--remove_taxa", ";".join(f"G{i}" for i in range(6)))
        self.assertEqual(p.returncode, 0, p.stderr)
        d = Path(self.dir)
        with (d / "t.decontam_pkg_decontaminated.csv").open() as f:
            self.assertEqual(list(csv.reader(f)), [["clade_name"] + [f"S{i}" for i in range(1, 13)]])
        for name in ["filter_report", "contaminants_final", "contaminants_by_batch"]:
            with (d / f"t.{name}.tsv").open() as f:
                self.assertEqual(f.readline().split("\t")[0], "taxon_id")

    def test_invalid_count_still_fails_even_when_all_would_be_filtered(self):
        d = Path(self.dir)
        text = (d / "otu.tsv").read_text().splitlines()
        fields = text[1].split("\t"); fields[1] = "-1"; text[1] = "\t".join(fields)
        (d / "otu.tsv").write_text("\n".join(text) + "\n")
        p = self.run_r("--min_abundance", "1000")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("Negative count values", p.stderr)
        self.assertFalse((d / "t.decontam_pkg_decontaminated.csv").exists())


class DependencyProbe(unittest.TestCase):
    def test_unavailable_dependency_skips_but_broken_R_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            executable = Path(tmp) / "Rscript"
            for code, expected in [(42, "False"), (0, "True"), (1, None)]:
                executable.write_text(f"#!/bin/sh\nexit {code}\n"); executable.chmod(0o700)
                p = subprocess.run([sys.executable, "-c", "import test_decontam_zero_total as z; print(z.decontam_available())"],
                                   cwd=ROOT / "tests", env=dict(os.environ, CMP_TEST_RSCRIPT=str(executable)),
                                   capture_output=True, text=True)
                if expected is None:
                    self.assertNotEqual(p.returncode, 0)
                    self.assertIn("R dependency probe failed", p.stderr)
                else:
                    self.assertEqual(p.returncode, 0, p.stderr)
                    self.assertEqual(p.stdout.strip(), expected)


if __name__ == "__main__":
    unittest.main()
