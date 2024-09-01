#!/usr/bin/bash

ic50=splitted_data/gdsc_cell_line_ic50_test_fraction_0.1_id_997_seed_42.csv
gep=data/gene_expression/gdsc-rnaseq_gene-expression.csv
smi=data/smiles/gdsc.smi
gdict=data/2128_genes.pkl
smi_lan=single_pytorch_model/smiles_language
param=single_pytorch_model/weights/best_mse_paccmann_v2.pt
res=results
hparam=single_pytorch_model/model_params.json

python examples/IC50/test_paccmann_ori.py $ic50 $gep $smi $gdict $smi_lan $param $res $hparam
