#!/bin/bash

# Load Conda
source ~/.bashrc

conda activate metaphlan4_env
echo "conda env activated"
echo "Bowtie2 version:"
bowtie2 --version | head -1

cd "${CMP_OUT_DIR:?set CMP_OUT_DIR to the directory that holds (or will hold) ref/}"/ref/metaphlan
echo "cd metaphlan dir"

echo "Decompressing .bz2 files..."
# bzip2 -dk mpa_vJun23_CHOCOPhlAnSGB_202307_SGB.fna.bz2
# bzip2 -dk mpa_vJun23_CHOCOPhlAnSGB_202307_VSG.fna.bz2

echo "Building Bowtie2 index with 32 threads..."
bowtie2-build --large-index --threads 32 \
  mpa_vJun23_CHOCOPhlAnSGB_202307_VSG.fna \
  mpa_vJun23_CHOCOPhlAnSGB_202307

echo "DONE - Bowtie2 index built successfully"
echo "Final files:"
ls -lh
