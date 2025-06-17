#!/usr/bin/bash
export PYTHONNOUSERSITE=1

net=None
omics=_data/SANGER_RNA_TPM.csv
out=_data/SANGER_RNA_Lin.pickle
python process_cell.py -net $net -omics $omics -out $out

omics=_data/CCLE_RNA_TPM.csv
out=_data/CCLE_RNA_Lin.pickle
python process_cell.py -net $net -omics $omics -out $out

omics=_data/GDSC_RNA_Array.csv
out=_data/GDSC_RNA_Lin.pickle
python process_cell.py -net $net -omics $omics -out $out
 
omics=_data/EXP_TCGA.csv
out=_data/TCGA_RNA_Lin.pickle
python process_cell.py -net $net -omics $omics -out $out

omics=_data/EXP_TCGA_ComBat.csv
out=_data/TCGA_RNA_Lin_CB.pickle
train=_data/TCGA_RNA_Lin.pickle
python process_cell.py -net $net -omics $omics -out $out -train $train

