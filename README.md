# GCNPath2024

![_GCNPath](https://github.com/user-attachments/assets/13d9a4da-4efa-4548-bec6-91459141d176)

GCNPath is a graph-based deep learning model designed for predicting anticancer drug response. This model leverages pathway-pathway association (PPA) graphs, which are compressed from STRING and RegNetwork, as well as a GSVA pathway correlation network. The training of GCNPath is conducted using transcriptome data from the SANGER Cell Model Passports.

The GCNPath2024 directory is originally the subdirectory within ```_IC50_Prediction/do_better```, where the benchmark tests are implemented. The ```_IC50_Prediction``` directory serves as the root directory of the GCNPath project, encompassing the benchmark tests, as well as the preprocessing steps for cell lines, drugs, and ln(IC<sub>50</sub>) data.

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
You can install the Conda environment with the following command:
```conda env create -f GCNPath.yaml```

If the command above does not work, please install the necessary packages manually. The success of the installation and model training may depend on the compatibility between PyTorch, PyTorch Geometric, CUDA toolkit, and your GPU and operating system. We trained and tested our model on Ubuntu 20.04.5 LTS using an NVIDIA GeForce RTX 3090.

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

# Implementation

## 1. Processing of cell data
```
bash process_cell.sh

1-1. Compress RNA data from gene- to pathway-level using GSVA
# Rscript process_cell_gsva.R \
#     data/cell_data/SANGER_RNA_TPM_Filt.csv \
#     data/path_data/c2.cp.biocarta.v2023.1.Hs.entrez.gmt \
#     processed/cell_data_biocarta/SANGER_RNA_GSVA.csv

1-2. Transform GSVA pathway data into graph form
# python process_cell.py \
#     -omics processed/cell_data_biocarta/SANGER_RNA_GSVA.csv \
#     -net data/net_data_biocarta/STR9_Reg_Corr_KNN5.csv \
#     -out processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle

1-2. Process external data with the scaler fitted by train data (parameter -train)
# python process_cell.py \
#     -omics processed/cell_data_biocarta/CCLE_RNA_GSVA_BROAD.csv \
#     -net data/net_data_biocarta/STR9_Reg_Corr_KNN5.csv \
#     -out processed/cell_data_biocarta/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle \
#     -train processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle
```

## 2. Processing of drug data
```
bash process_drug.sh
# python process_drug.py \
#    -smi data/drug_data/SMILES_GDSC.csv \
#    -col_names Drug_CID -col_smi SMILES_CAN \
#    -out_dir processed/drug_data/GDSC_Drug_Custom.pickle \
```

## 3. Model Training

### 3-1. Model Training in Various Test Scenarios
Model training in outer cross-validation across different test scenarios is handled by the ```train.sh```, which sequentially executes ```train_write.sh``` and ```train.py```. The file ```train.sh``` contains the list of input file paths and hyperparameters. In a meanwhile, ```train_write.sh``` contains the resource management parameters of CPU, RAM and GPU via SLURM. This script generates new bash files in the ```exe``` folder (e.g. ```GCN0_N0_RGCN.sh```), incorporating all input file paths and hyperparameters, which are then used to execute the ```train.py``` script. If you utilize SLURM with ```use_slurm``` as 1 within  ```train.sh```, all log files will be created in the ```out``` folder.

The columns for a cell line, drug, and ln(IC<sub>50</sub>) can be specified using the ```-col_cell```, ```-col_drug``` and ```-col_ic50```, respectively. The train fold in cross-validation is corresponding to ```-nth```, whose range is [0, 24] in strict-blind tests or [0, 9] in the rest ones. ```train.sh``` takes the following two parameters.  
IC<sub>50</sub> data : 0 [GDSC1+2], 1 [GDSC1], 2 [GDSC2]  
Test type : 0 [Normal], 1 [Cell-Blind], 2 [Drug-Blind], 3 [Strict-Blind]

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
Training a model using the entire GDSC1+2 dataset without data splitting is managed by the ```retrain_total.sh```, which sequentially executes ```train_write.sh``` and ```retrain_total.py```. The overall process is similar to the one described in **Section 3-1**. 

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
Model testing is conducted using test bash scripts (e.g. ```test_ccle.sh```, ```test_tcga.sh```, ```test_chembl.sh```), which sequentially executes ```test_write.sh``` and ```test.py```. The process is similar to that described in **Section 3-1**. If you want to only output predicted response values without calculating performance metrics, set the parameter ```-col_ic50 0```.

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
