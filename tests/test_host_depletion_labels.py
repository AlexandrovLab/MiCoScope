#!/usr/bin/env python3
"""The --t2t_phix_db index is T2T-CHM13 + PhiX since 24 Nov 2026; docs and build scripts must agree.

Before: human-GCA-phix-db.mmi held T2T-CHM13v2.0 alone (25 sequences) because the build script never
downloaded phiX174, and these tests stopped any file from claiming PhiX removal (port review #19).
Now: human-T2T-phiX-db.mmi holds the 25 T2T sequences plus NC_001422.1 (Ludmil review #47/#37:
PhiX spike-in was most of the "synthetic construct" signal in CRC50). These tests stop the docs from
still saying "no phiX", and the shipped build script from silently dropping PhiX again.
Code and docs only -- the index is not read here.
"""
import importlib.util
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_cohort_report_labels", ROOT / "scripts" / "build_cohort_report.py")
CR = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CR)

# Places allowed to mention the old T2T-only index: the ones that explain the change.
EXPLAINERS = {"docs/CHANGELOG.md", "docs/host_depletion_references.md", "docs/outputs.md",
              "conf/databases.tsv", "tests/test_host_depletion_labels.py", "tests/test_database_versions.py"}
STALE = re.compile(r"contains\s+\*{0,2}no\s+phi\s*X|no\s+phi\s*X\s+despite", re.I)


class Labels(unittest.TestCase):
    def test_cohort_report_stage_label(self):
        labels = {key: label for key, label, _ in CR.STAGES}
        self.assertEqual(labels["pct_removed_t2t_phix"], "Removed: T2T-CHM13 + PhiX")

    def test_no_file_still_says_the_index_lacks_phix(self):
        offenders = []
        for path in sorted(ROOT.rglob("*")):
            rel = path.relative_to(ROOT).as_posix()
            if (not path.is_file() or rel.startswith(".git/") or rel in EXPLAINERS
                    or path.suffix not in (".md", ".py", ".nf", ".config", ".html", ".R", ".r", ".sh")):
                continue
            if STALE.search(path.read_text(errors="replace")):
                offenders.append(rel)
        self.assertEqual(offenders, [])

    def test_reference_doc_and_links(self):
        doc = (ROOT / "docs" / "host_depletion_references.md").read_text()
        self.assertIn("human-T2T-phiX-db.mmi", doc)
        self.assertIn("26 sequences", doc)
        self.assertIn("NC_001422.1", doc)
        for path in ("README.md", "docs/parameters.md", "docs/outputs.md"):
            with self.subTest(path=path):
                self.assertIn("host_depletion_references.md", (ROOT / path).read_text())

    def test_site_config_uses_the_phix_index(self):
        self.assertIn("human-T2T-phiX-db.mmi", (ROOT / "conf" / "tscc.config").read_text())

    def test_names_unchanged(self):
        """The parameter and the QC columns keep their names (now accurate)."""
        self.assertIn("params.t2t_phix_db", (ROOT / "main.nf").read_text())
        qc = (ROOT / "scripts" / "build_qc_summary.py").read_text()
        self.assertIn('"reads_after_t2t_phix"', qc)
        self.assertIn('"pct_removed_t2t_phix"', qc)

    def test_the_build_script_downloads_phix_and_stops_on_failure(self):
        """The code-level cause of the old T2T-only index, fixed: PhiX is downloaded, and a failed
        download or concatenation stops the script instead of silently indexing T2T alone."""
        script = (ROOT / "scripts" / "download_references.sh").read_text()
        code = "\n".join(l for l in script.splitlines() if not l.lstrip().startswith("#"))
        self.assertRegex(code, r"(wget|curl)[^\n]*NC_001422")
        self.assertRegex(code, r"\[ -s ref/NC_001422\.fna \] \|\| \{[^}]*exit 1")
        self.assertIn("ref/NC_001422.fna > ref/human-T2T-phiX.fna", code)
        index = (ROOT / "scripts" / "create_minimap2_indexes.sh").read_text()
        self.assertIn("-d ref/human-T2T-phiX-db.mmi ref/human-T2T-phiX.fna", index)


if __name__ == "__main__":
    unittest.main()
