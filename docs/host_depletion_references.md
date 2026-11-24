# What the host-depletion indexes contain

`--hg38_db` and `--t2t_phix_db` are the two minimap2 indexes every run depletes against before
anything is classified, so every MiCoScope count is what these two did not remove. On TSCC both are
the lab's shared files (paths in `conf/tscc.config`), the same two
files pksProfiler uses (except that MiCoScope's T2T pass also holds PhiX since 24 Nov 2026; see below). Both are rows in [`conf/databases.tsv`](../conf/databases.tsv), and every run
records whether they are still the frozen files (*DATABASE VERSIONS* in `RUN_REPORT.txt`).

The contents below were read out of the index files themselves (their header and sequence dictionary)
by pksProfiler on 2026-10-22 (`scripts/describe_minimap2_index.py`, documented in pksProfiler's
`docs/host_depletion_references.md`). They were not re-read for this page; `conf/databases.tsv`
carries pksProfiler's dictionary sha256 for each index as its key, so a run against a different index
reports `CHANGED: key`.

## `human-GRC-db.mmi` — `--hg38_db`

GRCh38.p14, the RefSeq full assembly (GCF_000001405.40): 705 sequences, 3,298,430,636 bp, built
`-x sr` (k=21, w=11). Primary chromosomes, mitochondrion, and the unlocalized, unplaced, ALT and patch
scaffolds. **Not** in it: hs38d1 decoys, chrEBV, HLA allele contigs, phiX.

## `human-T2T-phiX-db.mmi` — `--t2t_phix_db`: T2T-CHM13 plus PhiX (since 24 Nov 2026)

T2T-CHM13v2.0 (GCA_009914755.4, 25 sequences, 3,117,292,070 bp: chr1–22, X, Y (HG002) and the
mitochondrion) **plus phiX174 (NC_001422.1, 5,386 bp)**: 26 sequences, 3,117,297,456 bp, indexed with
`minimap2 -x sr` (2.31). The FASTAs come from `scripts/download_references.sh` and the index from
`scripts/create_minimap2_indexes.sh`. The TSCC build (24 Nov 2026) also
checked the T2T FASTA's md5 against NCBI's checksum file, PhiX's length (5,386 bp), that the combined FASTA
has exactly 26 sequences and that length, and that every 150 bp tile of PhiX maps to NC_001422.1 (105 of 105). sha256 values are in `human-T2T-phiX-db.sha256`
beside the index; `conf/databases.tsv` freezes it.

**Why it was rebuilt.** The earlier file, `human-GCA-phix-db.mmi` (2025-01-19), held the 25 T2T sequences
only, despite its name: `scripts/download_references.sh` (copied from cguccione/human_host_filtration)
runs `cat ref/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna ref/NC_001422.fna > ref/human-GCA-phix.fna` but
never downloads `NC_001422.fna`, and without `set -e` the missing file only printed an error. PhiX
spike-in reads therefore reached KrakenUniq, whose database files UniVec's PhiX entry under **synthetic
construct** (taxid 32630) — not under viruses, as this page previously assumed. In CRC50 (Ludmil review
#47, 24 Nov 2026) PhiX was 91% of the reads matching UniVec, concentrated in one sequencing centre where
synthetic construct was ~98.5% of all classified reads; the rest were index and PCR primer fragments.

**What changes with the new index.** PhiX spike-in is removed at host depletion, before classification,
so it no longer counts in Bracken totals, library sizes or the synthetic-construct share. Read
`reads_after_t2t_phix` as "reads surviving the T2T + PhiX pass". Results change for libraries that carried
a PhiX spike-in; the CRC50 before/after regression (#8) measures it. Primer and adapter fragments are not
PhiX and are still classified as synthetic construct, which the microbial-scope gate (#14) and the
analysis tables (#21, #37) already exclude.

The old `human-GCA-phix-db.mmi` is kept for reproducing earlier runs (pksProfiler also still reads it).

## Re-checking

Identity, cheap (a few kB of each index; what `conf/databases.tsv` freezes as `mmi:dictionary`): every
run does it. Integrity, expensive (reads the whole 9 GB file; run deliberately, not per run):

```
sha256sum <db_dir>/human-GRC-db.mmi
sha256sum <db_dir>/human-T2T-phiX-db.mmi   # compare with human-T2T-phiX-db.sha256
```

---
Back to the [README](../README.md).
