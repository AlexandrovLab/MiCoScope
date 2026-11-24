#!/bin/bash


# Load Conda
source ~/.bashrc

# Utility script, not part of the pipeline run. Set TOOLS_ENV to an env that has the
# tools below; `microbiome` is the historical name and does not exist for every user.
conda activate "${TOOLS_ENV:-microbiome}" || { echo "ERROR: conda env '${TOOLS_ENV:-microbiome}' not found; rerun with TOOLS_ENV=<env with the tools this script needs>" >&2; exit 1; }
echo "conda env activated"

cd "${CMP_OUT_DIR:?set CMP_OUT_DIR to the directory that holds (or will hold) ref/}"
echo "cd ref dir"

bash ../scripts/download_references.sh

echo "DONE"
