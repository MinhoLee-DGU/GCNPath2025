# GCNPath2024
GCNPath is a deep learning model of anticancer drug response prediction. This model utilizes pathway-pathway association (PPA) graphs compressed from STRING and RegNetwork and a GSVA pathway correlation network. This model is trained with cell-line transcriptome from SANGER Cell Passports.

The GCNPath2024 directory is originally the subdirectory of ```_IC50_Prediction/do_better```, where the benchmark tests are implemented. ```_IC50_Prediction``` is the root directory of GCNPath project including the benchmark tests and preprocessing of cell-lines, drugs and ln(IC50) data.


# Quick start
```
conda env create -f GCNPath.yaml
conda activate GCNPath
bash process_cell.sh
bash process_drug.sh
bash train.sh 0 0
bash test_ccle.sh
```

# Requirement
You can install the environment via anaconda as follows.  
```conda env create -f GCNPath.yaml```

If the command above does not work, please install the necessary packages manually. The sucess of installation and training might depend on the compatibility between PyTorch, Pytorch Geometric, Cudatoolkit and your GPU and OS system. We trained and tested our model in Ubuntu 20.04.5 LTS using NVIDIA GeForce RTX3090.

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
#     data/cell_data/SANGER_RNA_TPM.csv \
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
Training a model in all fold is implemented by ```train.sh```, which sequentially executes ```train_write.sh``` and ```train.py```. The file ```train.sh``` contains the list of input files and hyper-parameters. In a meanwhile, ```train_write.sh``` contains the resource management parameters of CPU, RAM and GPU via SLURM. This file writes all inputs and parameters into new bash files in ```exe``` folder, which eventually implement the ```train.py```. All log files will be created in the ```out``` folder. 

The column of cell-line, drug, ln(IC50) can be designated with ```-col_cell```, ```-col_drug``` and ```-col_ic50```, respectively. The train fold is corresponding to ```-nth```, whose range is [0, 24] in strict-blind tests or [0, 9] in the rest ones. ```train.sh``` takes the following two parameters.  
IC50 data : 0 [GDSC1+2], 1 [GDSC1], 2 [GDSC2]  
Test type : 0 [Normal], 1 [Cell-Blind], 2 [Drug-Blind], 3 [Strict-Blind]

```
bash train.sh 0 0
# python train.py \
#    -cell processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle \
#    -drug processed/drug_data/GDSC_Drug_Custom.pickle \
#    -ic50 data/ic50_data/IC50_GDSC.txt \
#    -out_dir results/IC50_GDSC/Normal/RGCN -nth 0 \
#    -col_cell Cell -col_drug Drug -col_ic50 IC50
```

## 4. Model Testing
Testing a model in all fold is implemented by test bash files (ex. ```test_ccle.sh```, ```test_tcga.sh```, ```test_chembl.sh```), which sequentially executes ```test_write.sh``` and ```test.py```. The description is the same as **3. Model Training**. If you want to only output predicted response values without performance, set the parameter ```-col_ic50 0```.

```
bash test_ccle.sh
# python test.py \
#    -cell processed/cell_data_biocarta/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle \
#    -drug processed/drug_data/CCLE_Drug_Custom.pickle \
#    -ic50 data/ic50_data/IC50_CCLE.txt \
#    -dir_param results/IC50_GDSC/Normal/RGCN/param_0.pt \
#    -dir_hparam results/IC50_GDSC/Normal/RGCN/hyper_param_0.pickle \
#    -out_file results/IC50_GDSC/Normal/RGCN/pred_ccle_0.csv \
#    -col_cell Cell_BROAD_ID -col_drug Drug_CID -col_ic50 LN_IC50
```
