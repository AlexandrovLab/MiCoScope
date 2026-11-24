#!/usr/bin/env python3
"""Convert hierarchical Bracken MPA to exact-ID rank tables and non-overlapping totals."""
import argparse
import csv
from count_tables import read_counts, write_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out_genus", required=True)
    ap.add_argument("--out_species", required=True)
    ap.add_argument("--out_library", required=True)
    a = ap.parse_args()
    samples, genus = read_counts(a.genus, "g")
    species_samples, species = read_counts(a.species, "s")
    if set(samples) != set(species_samples):
        raise ValueError("Bracken genus/species sample IDs differ")
    # Exclusions applied to every analysis table, on both entry routes: viruses (the existing
    # consensus rule), and human and synthetic construct.
    # Human and vector reads are not microbes of the sample, so they are not analysis rows; they
    # stay visible in the QC summary's kraken_/bracken_ columns, and the raw per-sample Bracken
    # reports remain published unfiltered. The library totals below are computed AFTER these
    # exclusions, so decontam's absence evidence is judged on microbial reads.
    def nonmicrobial(clade: str) -> bool:
        tail = clade.rsplit("|", 1)[-1]
        name = tail.split("__", 1)[-1].replace("_", " ").strip().lower()
        return name in {"homo sapiens", "homo", "synthetic construct"}

    genus = {t: r for t, r in genus.items() if "Viruses" not in t and not nonmicrobial(t)}
    species = {t: r for t, r in species.items() if "Viruses" not in t and not nonmicrobial(t)}
    write_counts(a.out_genus, samples, genus)
    write_counts(a.out_species, samples, species)
    with open(a.out_library, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "bracken_genus_total", "bracken_species_total"])
        for s in samples:
            writer.writerow([s, sum(r[s] for r in genus.values()), sum(r[s] for r in species.values())])


if __name__ == "__main__":
    main()
