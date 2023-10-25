#!/bin/bash

# activate conda: https://github.com/conda/conda/issues/7126
. /fh/fast/bloom_j/software/miniconda3/etc/profile.d/conda.sh
conda activate ZIKV_DMS_NS3_EvansLab

snakemake -j 72 --keep-incomplete
