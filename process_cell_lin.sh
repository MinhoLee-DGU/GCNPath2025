#!/usr/bin/bash
export PYTHONNOUSERSITE=1

net=None
omics=processed/cell_data_biocarta/SANGER_RNA_GSVA.csv
out=processed/cell_data_biocarta/SANGER_RNA_Lin.pickle

python process_cell.py -net $net -omics $omics -out $out
