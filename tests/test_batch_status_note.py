#!/usr/bin/env python3
"""Batch sample ledger notes reach the cohort QC summary.

An included sample keeps the plain status 'included' when its reason is 'eligible', and shows the
recorded note otherwise (mixed sample types without sample type among the covariates). Excluded
decontamination controls keep their reason.
"""

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("build_qc_summary", ROOT / "scripts" / "build_qc_summary.py")
QC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(QC)

NOTE = "sample types mixed (Primary Tumor, Solid Tissue Normal); sample type is not a protected covariate"


class BatchStatusNote(unittest.TestCase):
    def ledger(self, rows):
        tmp = Path(tempfile.mkdtemp())
        path = tmp / "species.batch_sample_status.tsv"
        path.write_text("sample\tstatus\treason\n" + "".join(f"{s}\t{st}\t{r}\n" for s, st, r in rows))
        return QC.parse_batch_sample_status(path)

    def test_eligible_stays_plain(self):
        self.assertEqual(self.ledger([("T1", "included", "eligible")]), {"T1": "included"})

    def test_note_is_shown(self):
        self.assertEqual(self.ledger([("T1", "included", NOTE)]), {"T1": f"included: {NOTE}"})

    def test_control_exclusion_keeps_reason(self):
        reason = "decontamination control (pass --include_controls to analyse)"
        self.assertEqual(self.ledger([("B1", "excluded", reason)]), {"B1": f"excluded: {reason}"})


if __name__ == "__main__":
    unittest.main()
