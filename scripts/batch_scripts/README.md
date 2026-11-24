# batch_scripts

## Safe to submit

- `main.sbatch` — launches the pipeline with `-profile tscc`
- `resume_main.sbatch` — the same, kept under the name older notes use

        cd <your CMPipeline checkout>
        mkdir -p logs
        sbatch scripts/batch_scripts/main.sbatch samples.csv

  Both take the sample sheet as `$1` and pass any further arguments to Nextflow. They
  run in `$SLURM_SUBMIT_DIR`, so they operate on the checkout you submit from.

## Run by hand

The six other `*.sbatch` files here are one-time setup jobs. They **have no `#SBATCH` directives**
(stripped in commit `ca7b846`), so Slurm would apply its 1 hour default walltime and kill a
download or index build partway; the right walltime and memory for each have not been measured.
They take their working directory from an environment variable and stop if it is unset:

| script | working directory |
|---|---|
| build_metaphlan_index.sbatch | `$CMP_OUT_DIR/ref/metaphlan` |
| create_minimap2_indexes.sbatch | `$CMP_OUT_DIR` |
| download_microbial_db.sbatch | `$CMP_OUT_DIR/ref` |
| download_references.sbatch | `$CMP_OUT_DIR` |
| unpack_microbial_dbs.sbatch | `$CMP_OUT_DIR/ref` |
| run_consensus_taxa.sbatch | `$CMP_RUN_DIR` (a run directory with `scripts/` and `RESULTS/`) |

Run them with explicit `sbatch` arguments, e.g.

    CMP_OUT_DIR=/path/to/out sbatch -A <account> -p <partition> -q <qos> -c 4 --mem 32G -t 24:00:00 \
        --chdir "$PWD" scripts/batch_scripts/download_references.sbatch

Directives added to these files must sit directly under the shebang: Slurm stops reading them
at the first executable line. Validate with `sbatch --test-only`.
