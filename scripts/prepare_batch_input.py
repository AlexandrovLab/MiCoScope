#!/usr/bin/env python3
"""Route a validated nonempty table to ConQuR, or report an explicit empty-input skip.

No corrected matrix, normalization or statistical result is fabricated for a skip.
Malformed counts/metadata remain errors, including when the count table has no taxa.
"""
import argparse
import csv
import shutil
from pathlib import Path

from count_tables import read_counts, write_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--rank", choices=("species", "genus"), required=True)
    ap.add_argument("--taxon-column", default="clade_name")
    ap.add_argument("--sample-column", required=True)
    ap.add_argument("--batch-column", required=True)
    ap.add_argument("--type-column", required=True)
    a = ap.parse_args()
    with open(a.input) as f:
        header = f.readline()
    delimiter = "\t" if "\t" in header else ","
    samples, taxa = read_counts(a.input, a.rank[0], delimiter, a.taxon_column)
    if taxa:
        # Preserve the exact nonempty input; correction retains its existing validation.
        shutil.copyfile(a.input, "ready_counts.tsv")
        return
    with open(a.metadata) as f:
        header = f.readline()
        f.seek(0)
        reader = csv.DictReader(f, delimiter="\t" if "\t" in header else ",")
        columns = reader.fieldnames or []
        required = {a.sample_column, a.batch_column, a.type_column}
        if len(set(columns)) != len(columns) or not required.issubset(columns):
            raise ValueError("Metadata requires unique columns including sample, batch and type")
        seen = set()
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed metadata row")
            if any(not row[k] or row[k] != row[k].strip() or row[k] in ("NA", "NaN") for k in required):
                raise ValueError("Metadata has missing or padded sample/batch/type values")
            sample = row[a.sample_column]
            if sample in seen:
                raise ValueError("Duplicate metadata sample: " + sample)
            seen.add(sample)
        if not set(samples).issubset(seen):
            raise ValueError("Count-table sample missing from metadata")
    reason = "no_taxa_in_input"
    with open(f"{a.rank}.batch_sample_status.tsv", "w", newline="") as f:
        out = csv.writer(f, delimiter="\t", lineterminator="\n")
        out.writerow(["sample", "status", "reason"])
        out.writerows((s, "skipped", reason) for s in samples)
    write_counts("empty_counts.csv", samples, {}, delimiter=",", taxon_column=a.taxon_column)
    Path("batch_skip_status.tsv").write_text(
        "stage\tstatus\treason\tn_samples\tn_taxa\n"
        f"batch_correction\tskipped\t{reason}\t{len(samples)}\t0\n"
        f"normalization\tskipped\t{reason}\t{len(samples)}\t0\n"
        f"permanova\tskipped\t{reason}\t{len(samples)}\t0\n")


if __name__ == "__main__":
    main()
