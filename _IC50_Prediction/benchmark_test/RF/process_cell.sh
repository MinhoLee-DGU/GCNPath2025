#!/usr/bin/bash
export PYTHONNOUSERSITE=1

# omics=_data/SANGER_RNA_TPM.csv
# out=_data/SANGER_RNA.pickle
# python process_cell.py -omics $omics -out $out
# 
# omics=_data/CCLE_RNA_TPM.csv
# out=_data/CCLE_RNA.pickle
# python process_cell.py -omics $omics -out $out
# 
# omics=_data/GDSC_RNA_Array.csv
# out=_data/GDSC_RNA.pickle
# python process_cell.py -omics $omics -out $out
 # 
# omics=_data/EXP_TCGA.csv
# out=_data/TCGA_RNA.pickle
# python process_cell.py -omics $omics -out $out
# 
# omics=_data/EXP_TCGA_ComBat.csv
# out=_data/TCGA_RNA_CB.pickle
# python process_cell.py -omics $omics -out $out

omics=_data/SANGER_RNA_TPM.csv
out=_data/SANGER_RNA_Noise.pickle
python process_cell.py -omics $omics -out $out
