#!/usr/bin/bash

ic50=$1

case $ic50 in 
	0) IC50=./_data/IC50_GDSC.csv ;;
	1) IC50=./_data/IC50_GDSC1.csv ;;
	2) IC50=./_data/IC50_GDSC2.csv ;;
esac

ic50=$IC50
gep=data/gene_expression/gdsc-rnaseq_gene-expression.csv
smi=_data/SMILES_GDSC.smi
gdict=data/2128_genes.pkl
smi_lan=single_pytorch_model/smiles_language
# param=single_pytorch_model/weights/best_mse_paccmann_v2.pt
param=single_pytorch_model/weights/best_pearson_paccmann_v2.pt
res=results
hparam=single_pytorch_model/model_params.json

source activate paccmann_predictor
python examples/IC50/test_paccmann.py $ic50 $gep $smi $gdict $smi_lan $param $res $hparam
