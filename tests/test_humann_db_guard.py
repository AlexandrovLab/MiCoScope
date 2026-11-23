"""HUMAnN / MetaPhlAn-database compatibility guard (review #15, 23 Nov 2026).

HUMAnN 3.9 refuses MetaPhlAn profiles not generated with a v3/vJun23-era database
(probed: a vJan26 profile fails with "not generated with the database version v3 or
vJun23"). main.nf must therefore warn and SKIP functional profiling when it is
requested together with a newer --metaphlan_index, instead of failing at the HUMAnN
stage hours into a run. These are wiring tests over main.nf's source, in the style of
the other module-script tests: the guard must exist, be applied to the resolved
skip_humann3 flag, accept the vJun23/v3 families, and skip (never error) otherwise.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MAIN = (ROOT / "main.nf").read_text()


def guard_body():
    m = re.search(r"def humannDbGuard\(boolean skip\) \{(.*?)\n\}", MAIN, re.S)
    if m is None:
        raise AssertionError("humannDbGuard(boolean skip) not found in main.nf")
    return m.group(1)


class HumannDbGuard(unittest.TestCase):
    def test_guard_wraps_the_resolved_skip_flag(self):
        self.assertIn(
            "humannDbGuard(boolParam('skip_humann3', params.skip_humann3))", MAIN,
            "skip_humann3 must be resolved through the guard, so every downstream "
            "consumer (gating, provenance, RUN_REPORT) sees the guarded value",
        )

    def test_guard_is_defined_outside_the_workflow_block(self):
        self.assertLess(
            MAIN.index("def humannDbGuard"), MAIN.index("\nworkflow {"),
            "the guard body must live outside the workflow block (v1-parser 65535-char limit)",
        )

    def test_vjun23_and_v3_families_pass_through(self):
        # Compatibility is defined once, in lib/Mpa.groovy, and the guard uses it.
        self.assertIn("Mpa.humannCompatible(", guard_body())
        mpa = (ROOT / "lib" / "Mpa.groovy").read_text()
        rule = re.search(r"static boolean humannCompatible\(String index\) \{(.*?)\n    \}", mpa, re.S)
        self.assertIsNotNone(rule, "Mpa.humannCompatible not found")
        self.assertIn("vJun23", rule.group(1))
        self.assertIn("mpa_v3", rule.group(1))

    def test_guard_reads_the_resolved_index(self):
        # With --metaphlan_index auto the guard must see the index Mpa resolved (vJun23 after the
        # HUMAnN switch), not the literal 'auto', or every HUMAnN run would be skipped.
        self.assertIn("Mpa.index(params)", guard_body())

    def test_incompatible_index_warns_and_skips_never_errors(self):
        body = guard_body()
        self.assertIn("log.warn", body, "must warn, visibly, at launch")
        self.assertIn("return true", body, "must skip (return true), letting the run continue")
        self.assertNotIn("error ", body, "must not fail the run; the user chooses by flag")
        self.assertIn("vJun23_CHOCOPhlAnSGB_202307", body,
                      "the warning must name the database that restores functional profiles")

    def test_workflow_block_stays_under_the_v1_parser_limit(self):
        m = re.search(r"\nworkflow\s*\{", MAIN)
        depth, end = 0, None
        for i in range(m.end() - 1, len(MAIN)):
            if MAIN[i] == "{":
                depth += 1
            elif MAIN[i] == "}":
                depth -= 1
                if depth == 0:
                    end = i
                    break
        self.assertLessEqual(end - (m.end() - 1) + 1, 65535,
                             "workflow block exceeds the Nextflow v1 parser limit")


if __name__ == "__main__":
    unittest.main()
