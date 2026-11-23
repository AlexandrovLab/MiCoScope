"""MetaPhlAn database / index / environment resolution, lib/Mpa.groovy (review #15, Nov 2026).

The default database is the newest the site configures (vJan26 on TSCC) with MetaPhlAn 4.2.
HUMAnN 3.9 reads only vJun23-era profiles, so --skip_humann3 false with --metaphlan_index auto
switches to --metaphlan_db_vjun23 and MetaPhlAn 4.1.1. Explicit settings win. --metaphlan_index
auto takes the index from the single mpa_*.pkl in the database directory, so callers that pass an
older database keep their old behaviour.

Each case runs the REAL lib/Mpa.groovy through Nextflow on placeholder database directories (empty
.pkl files only). Skipped when Nextflow is not available (set NEXTFLOW to point at it).
"""
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NF_ENV = Path("/tscc/projects/ps-lalexandrov/amabbasi/anaconda3/envs/nf-env")


def find_nextflow():
    cand = os.environ.get("NEXTFLOW") or shutil.which("nextflow")
    if not cand and (NF_ENV / "bin" / "nextflow").exists():
        cand = str(NF_ENV / "bin" / "nextflow")
    return cand


NEXTFLOW = find_nextflow()

HARNESS = r'''
params.metaphlan_db = null
params.metaphlan_db_vjun23 = null
params.metaphlan_index = 'auto'
params.metaphlan4_env = 'auto'
params.metaphlan41_env = 'ENV41'
params.metaphlan42_env = 'ENV42'
params.skip_humann3 = true
workflow {
    def err = Mpa.check(params)
    def db = Mpa.dbDir(params)
    println "RESULT|${Mpa.switched(params)}|${db ? new File(db).name : 'null'}|${Mpa.index(params)}|${Mpa.env(params)}|${err ? 'ERR:' + err : 'ok'}"
}
'''


@unittest.skipUnless(NEXTFLOW, "nextflow not available")
class MpaResolution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="mpa_res_"))
        (cls.tmp / "lib").mkdir()
        shutil.copy(ROOT / "lib" / "Mpa.groovy", cls.tmp / "lib" / "Mpa.groovy")
        (cls.tmp / "main.nf").write_text(HARNESS)
        for d, pkls in {"db_new": ["mpa_vJan26_CHOCOPhlAnSGB_202605.pkl"],
                        "db_old": ["mpa_vJun23_CHOCOPhlAnSGB_202307.pkl"],
                        "db_two": ["mpa_a.pkl", "mpa_b.pkl"]}.items():
            (cls.tmp / d).mkdir()
            for p in pkls:
                (cls.tmp / d / p).write_text("")
        cls.env = dict(os.environ, NXF_ANSI_LOG="false", NXF_VER=os.environ.get("NXF_VER", "25.10.4"))
        if (NF_ENV / "lib" / "jvm").exists():
            cls.env["JAVA_HOME"] = str(NF_ENV / "lib" / "jvm")
            cls.env["PATH"] = f"{NF_ENV / 'lib' / 'jvm' / 'bin'}:{cls.env.get('PATH', '')}"
        # One warm-up launch: the first Nextflow run in a fresh directory initialises its state, and
        # the first test case once failed there without reproducing in five later runs (23 Nov 2026).
        subprocess.run([NEXTFLOW, "-q", "run", "main.nf", "--metaphlan_db", str(cls.tmp / "db_old")],
                       cwd=cls.tmp, env=cls.env, capture_output=True, text=True, timeout=600)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def resolve(self, *args):
        flags = []
        for a in args:
            flags.extend(a.replace("{T}", str(self.tmp)).split())
        r = subprocess.run([NEXTFLOW, "-q", "run", "main.nf", *flags], cwd=self.tmp, env=self.env,
                           capture_output=True, text=True, timeout=600)
        line = next((l for l in r.stdout.splitlines() if l.startswith("RESULT|")), None)
        self.assertIsNotNone(line, f"no RESULT line; stderr:\n{r.stderr[-2000:]}")
        _, switched, db, index, env, status = line.split("|", 5)
        return dict(switched=switched, db=db, index=index, env=env, status=status)

    def test_default_uses_the_newest_database_and_4_2(self):
        r = self.resolve("--metaphlan_db {T}/db_new --metaphlan_db_vjun23 {T}/db_old")
        self.assertEqual((r["switched"], r["db"], r["index"], r["env"], r["status"]),
                         ("false", "db_new", "mpa_vJan26_CHOCOPhlAnSGB_202605", "ENV42", "ok"))

    def test_humann_switches_to_vjun23_and_4_1(self):
        r = self.resolve("--metaphlan_db {T}/db_new --metaphlan_db_vjun23 {T}/db_old --skip_humann3 false")
        self.assertEqual((r["switched"], r["db"], r["index"], r["env"], r["status"]),
                         ("true", "db_old", "mpa_vJun23_CHOCOPhlAnSGB_202307", "ENV41", "ok"))

    def test_humann_with_an_old_database_given_directly_needs_no_switch(self):
        r = self.resolve("--metaphlan_db {T}/db_old --skip_humann3 false")
        self.assertEqual((r["switched"], r["db"], r["env"], r["status"]), ("false", "db_old", "ENV41", "ok"))

    def test_humann_without_a_vjun23_database_stops_at_launch_naming_the_parameter(self):
        r = self.resolve("--metaphlan_db {T}/db_new --skip_humann3 false")
        self.assertTrue(r["status"].startswith("ERR:"), r)
        self.assertIn("--metaphlan_db_vjun23", r["status"])

    def test_explicit_index_wins_over_the_humann_switch(self):
        r = self.resolve("--metaphlan_db {T}/db_new --metaphlan_db_vjun23 {T}/db_old --skip_humann3 false "
                         "--metaphlan_index mpa_vJan26_CHOCOPhlAnSGB_202605")
        self.assertEqual((r["switched"], r["db"], r["index"]), ("false", "db_new", "mpa_vJan26_CHOCOPhlAnSGB_202605"))

    def test_explicit_env_wins(self):
        r = self.resolve("--metaphlan_db {T}/db_new --metaphlan4_env /my/env.txt")
        self.assertEqual(r["env"], "/my/env.txt")

    def test_ambiguous_database_directory_stops_and_asks_for_an_index(self):
        r = self.resolve("--metaphlan_db {T}/db_two")
        self.assertTrue(r["status"].startswith("ERR:"), r)
        self.assertIn("--metaphlan_index", r["status"])

    def test_legacy_callers_with_an_old_database_are_unchanged(self):
        r = self.resolve("--metaphlan_db {T}/db_old")
        self.assertEqual((r["switched"], r["index"], r["env"], r["status"]),
                         ("false", "mpa_vJun23_CHOCOPhlAnSGB_202307", "ENV41", "ok"))


if __name__ == "__main__":
    unittest.main()
