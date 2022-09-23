#!/usr/bin/bash

nextflow run main.nf \
    --outdir data/results \
    -resume \
    -profile icbi_liver \
    -w /data/scratch/sturm/projects/2022/schneeberger-liver/work
