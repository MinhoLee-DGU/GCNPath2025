#!/usr/bin/bash
python examples/IC50/test_paccmann_ori.py \
    splitted_data/gdsc_example_322.csv \
    data/gene_expression/gdsc-rnaseq_gene-expression.csv \
    data/smiles/gdsc.smi \
    data/2128_genes.pkl \
    single_pytorch_model/smiles_language \
    single_pytorch_model/weights/best_mse_paccmann_v2.pt \
    results_bs128 \
    single_pytorch_model/model_params.json
