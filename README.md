# GCNPath2025

![_GCNPath](https://github.com/user-attachments/assets/13d9a4da-4efa-4548-bec6-91459141d176)

GCNPath is a graph-based deep learning model designed for predicting anticancer drug responses. The model utilizes pathway crosstalk network (PCN) graphs, which are compressed from STRING and RegNetwork, along with a GSVA pathway correlation network. GCNPath is trained using transcriptome data from the SANGER Cell Model Passports.

The GCNPath2025 directory was originally a subdirectory within ```_IC50_Prediction/do_better```, where the benchmark tests are implemented. The ```_IC50_Prediction``` directory serves as the root directory for the GCNPath project, containing benchmark tests as well as preprocessing steps for cell lines, drugs, and ln(IC<sub>50</sub>) data.

# Quick start
```
conda env create -f GCNPath.yaml
conda activate GCNPath
bash process_cell.sh
bash process_drug.sh
# bash train.sh 0 0
bash retrain_total.sh
bash test_ccle.sh
```

# Requirement
To install the Conda environment, run:
```
conda env create -f GCNPath.yaml
```

If the above command doesn't work, manually install the required packages. Installation and model training success depend on the compatibility between PyTorch, PyTorch Geometric, the CUDA toolkit, your GPU, and your operating system. We trained and tested the model on Ubuntu 20.04.5 LTS using an NVIDIA GeForce RTX 3090.
Required packages:
```
* python (3.8.18)
* numpy (1.24.4)
* pandas (2.0.3)
* scikit-learn (1.2.2)
* pytorch (1.11.0)
* torchaudio (0.11.0)
* torchvision (0.12.0)
* cudatoolkit (11.3.1)
* pyg (2.1.0)
* pytorch-cluster (1.6.0)
* pytorch-scatter (2.0.9)
* pytorch-sparse (0.6.15)
* rdkit (2022.09.5)
* dgl (1.1.0)
* dgl-life (0.2.9)
* libstdcxx-ng (13.2.0)
* r-base (4.2.3)
* r-matrixstats (1.1.0)
* bioconductor-gsva (1.46.0)
```

# Implementation

## 1. Cell Data Processing
```
bash process_cell.sh

1-1. Compress RNA data from gene-level to pathway-level using GSVA.
# Rscript process_cell_gsva.R \
#     data/cell_data/SANGER_RNA_TPM_Filt.csv \
#     data/path_data/c2.cp.biocarta.v2023.1.Hs.entrez.gmt \
#     processed/cell_data_biocarta/SANGER_RNA_GSVA.csv

1-2. Transform GSVA pathway data into graph format.
# python process_cell.py \
#     -omics processed/cell_data_biocarta/SANGER_RNA_GSVA.csv \
#     -net data/net_data_biocarta/STR9_Reg_Corr_KNN5.csv \
#     -out processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle

1-2. Process external data with the scaler fitted to training data (using the -train parameter).
# python process_cell.py \
#     -omics processed/cell_data_biocarta/CCLE_RNA_GSVA_BROAD.csv \
#     -net data/net_data_biocarta/STR9_Reg_Corr_KNN5.csv \
#     -out processed/cell_data_biocarta/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle \
#     -train processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle
```

## 2. Drug Data Processing
```
bash process_drug.sh
# python process_drug.py \
#    -smi data/drug_data/SMILES_GDSC.csv \
#    -col_names Drug_CID -col_smi SMILES_CAN \
#    -out_dir processed/drug_data/GDSC_Drug_Custom.pickle \
```

## 3. Model Training

### 3-1. Model Training in Various Test Scenarios
Model training in outer cross-validation across different test scenarios is handled by the ```train.sh``` script, which sequentially executes ```train_write.sh``` and ```train.py```. The ```train.sh``` file contains a list of input file paths and hyperparameters. Meanwhile, ```train_write.sh``` contains resource management parameters for CPU, RAM, and GPU via SLURM. This script generates new bash files in the ```exe``` folder (e.g., ```GCN0_N0_RGCN.sh```), incorporating all input file paths and hyperparameters, which are then used to execute the ```train.py``` script. If SLURM is used (with ```use_slurm``` set to 1 within ```train.sh```), log files will be created in the ```out``` folder.

The columns for cell lines, drugs, and ln(IC</sub>50</sub>) can be specified using ```-col_cell```, ```-col_drug```, and ```-col_ic50```, respectively. The training fold in cross-validation corresponds to ```-nth```, with a range of [0, 24] for strict-blind tests or [0, 9] for others. The ```train.sh``` script takes the following parameters:
```
* IC</sub>50</sub> data : 
** 0 [GDSC1+2]
** 1 [GDSC1]
** 2 [GDSC2]
* Test type : 
** 0 [Unblinded]
** 1 [Cell-Blind]
** 2 [Drug-Blind]
** 3 [Strict-Blind]
```

You can set the random seed for initializing model parameter weights using the ```-seed_model (default 2021)```. Note that the seed is used to assess the stability of model performance, rather than to reproduce the exact same prediction results. This is due to non-deterministic operations within PyTorch Geometric modules, such as ```torch_scatter```or when training models quickly using multiple workers for data loading with the ```-cpu```.

```
bash train.sh 0 0
# python train.py \
#    -cell processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle \
#    -drug processed/drug_data/GDSC_Drug_Custom.pickle \
#    -ic50 data/ic50_data/IC50_GDSC.txt \
#    -out_dir results/IC50_GDSC/Normal/RGCN -nth 0 \
#    -col_cell Cell -col_drug Drug -col_ic50 LN_IC50 -cpu 4
```

### 3-2. Model Training with Whole Dataset without Splitting Data
Training a model using the entire GDSC1+2 dataset without splitting data is managed by the ```retrain_total.sh``` script, which sequentially executes ```train_write.sh``` and ```retrain_total.py```. The overall process is similar to that described in **Section 3-1**. Since the GDSC1+2 dataset is not split into training, validation, and test sets, early stopping is not applied when training models for 100 epochs by default.

```
bash retrain_total.sh
# python retrain_total.py \
#    -cell processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle \
#    -drug processed/drug_data/GDSC_Drug_Custom.pickle \
#    -ic50 data/ic50_data/IC50_GDSC.txt \
#    -out_dir results/IC50_GDSC/Normal/RGCN \
#    -col_cell Cell -col_drug Drug -col_ic50 LN_IC50 -cpu 4 -seed_model 2021
```

## 4. Model Testing
Model testing is performed using test bash scripts (e.g., ```test_ccle.sh```, ```test_tcga.sh```, ```test_chembl.sh```), which sequentially execute ```test_write.sh``` and ```test.py```. The process is similar to **Section 3-1**. To output only predicted response values without calculating performance metrics, set the parameter ```-col_ic50``` to 0.

```
bash test_ccle.sh
# python test.py \
#    -cell processed/cell_data_biocarta/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle \
#    -drug processed/drug_data/CCLE_Drug_Custom.pickle \
#    -ic50 data/ic50_data/IC50_CCLE.txt \
#    -dir_param results/IC50_GDSC/Normal/RGCN/param_retrain_seed2021.pt \
#    -dir_hparam results/IC50_GDSC/Normal/RGCN/hyper_param_retrain_seed2021.pickle \
#    -out_file results/IC50_GDSC/Normal/RGCN/pred_ccle_seed2021.csv \
#    -col_cell Cell_BROAD_ID -col_drug Drug_CID -col_ic50 LN_IC50 -cpu 4
```
