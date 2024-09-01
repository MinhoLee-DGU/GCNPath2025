#!/usr/bin/bash

train_sensitivity_filepath="data/drug_sensitivity/gdsc-cell-line-name_drug-sensitivity_ex.csv"
test_sensitivity_filepath="data/drug_sensitivity/gdsc-cell-line-name_drug-sensitivity_ex.csv"
gep_filepath="data/gene_expression/gdsc-rnaseq_gene-expression.csv"
smi_filepath="data/smiles/gdsc.smi"
gene_filepath="data/2128_genes.pkl"
smiles_language_filepath="single_pytorch_model/smiles_language"
model_path="train_example"
params_filepath="single_pytorch_model/model_params.json"
training_name="train_example"

python examples/IC50/train_paccmann_ori.py \
    $train_sensitivity_filepath \
    $test_sensitivity_filepath \
    $gep_filepath \
    $smi_filepath \
    $gene_filepath \
    $smiles_language_filepath \
    $model_path \
    $params_filepath \
    $training_name
