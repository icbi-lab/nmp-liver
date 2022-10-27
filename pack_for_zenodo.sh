#!/bin/bash

# singularity containers needed to run all processes in the workflow
# 7za a data/zenodo/containers.zip \
#   containers/deseq2.sif \
#   containers/2022-schneeberger-liver-scanpy_2022-10-24.sif

# input data required to run all workflows from scratch
# 7za a data/zenodo/input_data.zip data/02_anndata

# final and intermediate results
7za a data/zenodo/results.zip data/results

