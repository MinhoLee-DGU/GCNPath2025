# GCNPath2025

![_GCNPath](https://github.com/user-attachments/assets/13d9a4da-4efa-4548-bec6-91459141d176)

GCNPath is a graph-based deep learning model designed for predicting anticancer drug responses. The model utilizes pathway crosstalk network (PCN) graphs, which are compressed from STRING and RegNetwork, along with a GSVA pathway correlation network. GCNPath is trained using transcriptome data from the SANGER Cell Model Passports.

The ```GCNPath2025``` directory was originally a subdirectory within ```_IC50_Prediction/benchmark_test```, where the benchmark tests are implemented. The ```_IC50_Prediction``` directory serves as the root directory for the GCNPath project, containing benchmark tests as well as preprocessing steps for cell lines, drugs, and ln(IC<sub>50</sub>) data.

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

# Directory Structure
## GCNPath2025
```
GCNPath2025
├── SMILESVec                       
│   └── source
│       └── process_drug_svec.sh    Process drug data into SMILESVec Fingerprint
├── _IC50_Prediction                Case study with cell lines (Colorectal, Breast, SCLC)
├── _case_study                     Case study with cell lines (Colorectal, Breast, SCLC)
├── _performance_tuning             Performances analysis in ablation tests
├── data                            Store daw data
│   ├── cell_data                    - Cell transcriptome data
│   ├── drug_data                    - Drug structure data in SMILES format
│   ├── ic50_data                    - IC50 or drug response data
│   ├── net_data_biocarta            - PCN graph structure data
│   └── path_data                    - Pathway annotation retrieved from MSigDB
├── exe                             Store executable bash files when running train.sh, retrain_total.sh and test_xxx.sh
├── model                           Model architectures
├── out                             Store SLURM error and log files when running train.sh, retrain_total.sh and test_xxx.sh
├── processed                       Store processed data
│   ├── cell_data_biocarta           - Store cell data processed with process_cell.sh
│   └── drug_data                    - Store drug data processed with process_drug.sh
├── results                         Store model weight parameters and prediction results
│   └── IC50_GDSC                    - Results when training with GDSC1+2 dataset and 4 test scenarios
│       ├── Unblinded
│       ├── Cell_Blind
│       ├── Drug_Blind
│       └── Strict_Blind
├── utils                           Utils
│   ├── functions.R                  - Utils for data analysis in _case_study and _performance_tuning
│   ├── utils.py                     - Utils for calculating performance, splitting data, and saving prediction results
│   ├── utils_gnn.py                 - Utils for graph neural networks
│   ├── utils_gsva.R                 - Utils for GSVA in process_cell_gsva.R
│   └── utils_model.py               - Utils for loading data, training and testing models, and model interpretability
├── GCNPath.yaml                    Conda environment file
├── process_cell.sh                 Process cell data into graph-based format (run process_cell_gsva.R and process_cell.py)
├── process_cell_lin.sh             Process cell data into non-graph-based format (run process_cell_gsva.R and process_cell.py)
├── process_cell_pert.sh            Process cell data with PCN graphs in perturbed topologies (run process_cell_gsva.R and process_cell.py)
├── process_cell_gsva.R             Process cell data from gene- into pathway-level with GSVA
├── process_cell.py                 Process cell data
├── process_drug.sh                 Process drug data into graph-based format (run process_drug.py)
├── process_drug_lin.sh             Process drug data into Morgan Fingerprint (run process_drug.py)
├── process_drug.py                 Process drug data
├── Readme.md                       README
├── train.sh                        Train model with GDSC datasets (25 or 10 outer cross-validation)
├── retrain_total.sh                Train model with GDSC datasets (no data split in train-valid-test)
├── train_write.sh                  Operated after train.sh & retrain_total.sh (contain SLURM configurations)
├── test_xxx.sh                     Test model with external datasets (e.g. ChEMBL, TCGA, SCLC)
└── test_write.sh                   Operated after test_xxx.sh (contain SLURM configurations)
```

## _IC50_Prediction
The ```GCNPath2025``` directory was originally a subdirectory within ```_IC50_Prediction/benchmark_test```
```
_IC50_Prediction
├── benchmark_test
│   ├── (SOTA model names)            SOTA models
│   ├── GCNPath                       GCNPath (You are here)
│   ├── _case_study_sclc              SCLC case study with Liu et al. (2024)
│   ├── _performance_chembl           Performance analysis with ChEMBL
│   ├── _performance_gdsc             Performance analysis with GDSC
│   ├── _performance_pretrain         Performance analysis with GDSC, using models pretrained by original model developers
│   ├── _performance_tcga             Performance analysis with TCGA
│   └── _performance_tcga_combat      Performance analysis with TCGA, after applying ComBat batch correction
├── processed_data
│   ├── cell_data                     Store cell data processed by "1-1_process_rna" and "4-2_tcga_response"
│   ├── drug_data                     Store drug data processed by "2-1_get_drug_cid"
│   ├── ic50_data                     Store ic50 data processed by "3-1_process_ic50" and "3-2_process_ic50_ccle"
│   └── net_data                      Store network data processed by "1-2_process_net" and "1-3_gsva_pcn"
├── project
│   ├── 1-1_process_rna               Process cell data from SANGER Cell Model Passports, CCLE and GDSC
│   ├── 1-2_process_net               Process network data (STRING, RegNetwork)
│   ├── 1-3_gsva_pcn                  Process PCN graphs
│   ├── 2-1_get_drug_cid              Process drug data from GDSC
│   ├── 3-1_process_ic50              Process IC50 data from GDSC
│   ├── 3-2_process_ic50_ccle         Process IC50 data from CCLE
│   ├── 4-1_external_chembl           Process IC50 data from ChEMBL
│   ├── 4-2_tcga_response             Process RNA-seq and drug response data from TCGA
│   ├── 4-3_SCLC_Liu_2024             Process RNA-seq and drug response data from Liu et al. (2024) for SCLC case study
│   ├── _prepare_(SOTA model names)   Process data to run SOTA models
│   ├── functions.R                   Utils for data analysis
│   └── utils_gsva.R                  Utils for GSVA in 1-3_gsva_pcn
└── raw_data                          Store raw data
```

# Requirement
To install the Conda environment, run:
```
conda env create -f GCNPath.yaml
```

If the above command doesn't work, manually install the required packages. Installation and model training success depend on the compatibility between PyTorch, PyTorch Geometric, the CUDA toolkit, your GPU and operating system. We trained and tested the model on Ubuntu 20.04.5 LTS using an NVIDIA GeForce RTX 3090. We are in the process of developing Docker and Singularity configurations to facilitate the setup of an environment based on Ubuntu 20.04 with Conda.

Required packages:
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

## 1. Processing Cell Data
Cell data are processed using ```process_cell.sh```, which sequentially executes ```process_cell_gsva.R``` and ```process_cell.py```. If PCN graphs are not provided (```-net None``` in ```process_cell.py```), the pathway score data are not formatted as graph[s], which can be implemented by the script ```process_cell_lin.sh```. We conducted an ablation test to evaluate the impact of PCN graph topology on prediction performance, using perturbed PCN graphs that preserve the degree of each pathway node, generated via ```process_cell_pert.sh```.

### process_cell_gsva.R
This script compresses RNA data from the gene level to the pathway level using GSVA. By default, genes not included in any pathway are filtered out.

* [```1st argument```] Gene-level transcriptome data in (input, CSV)
* [```2nd argument```] Pathway list (input, GMT)
* [```3rd argument```] Pathway score data. The file is structured with ```cell × pathway``` in row × column (output, CSV)
* [```4th argument```] Whether ```1st argument``` is structured with ```cell × gene``` in row × column (Default: TRUE)

### process_cell.py
This script standardizes pathway score data using ```RobustScaler``` and formats them into PCN graph[s]. When processing external cell data not used during training, pathway score data utilized in training stage must be provided to apply the standardization scaler for transforming these external data (```-train``` parameter).

* [```-omics```] Pathway score data processed in the previous step (input, CSV)
* [```-net```] PCN graphs containing at least two columns (```Pathway1```, ```Pathway2```, [```Edge_Type```]) (input, CSV)
* [```-out```] Pathway score data formatted as PCN graph[s] (output, Pickle)
* [```-train```] Output file containing the standardization scaler (optional input, Pickle)

```
bash process_cell.sh

1-1. Compress RNA data from gene-level to pathway-level using GSVA.
# Rscript process_cell_gsva.R \
#     data/cell_data/SANGER_RNA_TPM_Filt.csv \
#     data/path_data/c2.cp.biocarta.v2023.1.Hs.entrez.gmt \
#     processed/cell_data_biocarta/SANGER_RNA_GSVA.csv

1-2. Transform GSVA pathway score data into graph format.
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

## 2. Processing Drug Data
Drug data are processed using ```process_drug.sh```, which executes ```process_drug.py``` to convert drug structures into 2D graphs. When graph featurization is deactivated (```-drug_feat 0```), drug data are processed in Morgan Fingerprints (256-bit, radius 2), which can be implemented by the script ```process_drug_lin.sh```. We also generated SMILESVec features via ```SMILESVec/source/process_drug_svec.sh```, which were also found to underperformed compared to graph-based features in ablation tests.

### process_drug.py
* [```-smi```] Drug structure data in SMILES format containing at least one columns (```-col_smi```, [```-col_name```]) (input, CSV)
* [```-out_dir```] Drug structure data in graph format (output, Pickle)
* [```-col_smi```] Column name of drug structures in ```-smi``` file (Default: ```Drug_CID```)
* [```-col_name```] Column name of drug IDs or names in ```-smi``` file. If not provided (```None```), the SMILES strings themselves are used as drug IDs. (optional, Default: ```SMILES_CAN```)
* [```-drug_feat```] Whether drug data are processed in graph format? (Default: ```3```)

```
bash process_drug.sh
# python process_drug.py \
#    -smi data/drug_data/SMILES_GDSC.csv \
#    -col_names Drug_CID -col_smi SMILES_CAN \
#    -out_dir processed/drug_data/GDSC_Drug_Graph.pickle \
```

## 3. Training Models

### 3-1. Training Models in Various Test Scenarios
Training models in outer cross-validation across different test scenarios is handled by the ```train.sh``` script, which sequentially executes ```train_write.sh``` and ```train.py```. The ```train.sh``` file contains a list of input file paths and hyperparameters. Meanwhile, ```train_write.sh``` contains resource management parameters for CPU, RAM, and GPU via SLURM. This script generates new bash files in the ```exe``` folder (e.g., ```GCN0_N0_RGCN.sh```), incorporating all input file paths and hyperparameters, which are then used to execute the ```train.py``` script. If SLURM is used (with ```use_slurm``` set to ```1``` within ```train.sh```), log files will be created in the ```out``` folder. The training fold in cross-validation corresponds to ```-nth```, with a range of [```0, 24```] for strict-blind tests or [```0, 9```] for others. The ```train.sh``` script takes the following parameters:

### train.sh
* [```1st argument```] IC<sub>50</sub> data from GDSC1+2, GDSC1 or GDSC2 (choose one of ```0-2```)
* [```2nd argument```] Test type of Unblinded, Cell-Blind, Drug-Blind, and Strict-Blind tests (choose one of ```0-3```)

In ```train.py```, the columns for cell lines, drugs, and ln(IC<sub>50</sub>) in IC<sub>50</sub> data can be specified using ```-col_cell```, ```-col_drug```, and ```-col_ic50```, respectively. You can set the random seed for initializing model parameter weights using the ```-seed_model (default 2021)```. Note that the seed is used to assess the stability of model performance, rather than to reproduce the exact same prediction results. This is due to non-deterministic operations within PyTorch Geometric modules, such as ```torch_scatter```or when training models quickly using multiple workers for data loading with the ```-cpu```.

### train.py
* [```-cell```] Pathway score data formatted as PCN graph[s] (output, Pickle)
* [```-drug```] Drug structure data in graph format (output, Pickle)
* [```-ic50```] IC50 data containing at least three columns (```-col_cell```, ```-col_drug```, ```-col_ic50```) (input, TXT or CSV)
* [```-out_dir```] Directory to save hyperparameter, weight parameter and prediction results (output, Directory)
* [```-nth```] Fold index for outer cross validation (choose one of ```0-24``` in Strict-Blind or ```0-9``` in the other test types)
* [```-col_cell```] Column name of cell in ```-ic50``` file (Default: ```Cell```)
* [```-col_drug```] Column name of drug in ```-ic50``` file (Default: ```Drug```)
* [```-col_ic50```] Column name of IC<sub>50</sub> in ```-ic50``` file (Default: ```LN_IC50```)
* [```-seed_model```] Seed number for model weight initialization (Default: ```2021```)
* [```-e```] Maximum epochs for training model (Default: ```300```)
* [```-p```] Patience of early stopping in training model (Default: ```10```)
* [```-cpu```] Number of workers for data loading (Default: ```4```).

```
bash train.sh 0 0
# python train.py \
#    -cell processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle \
#    -drug processed/drug_data/GDSC_Drug_Graph.pickle \
#    -ic50 data/ic50_data/IC50_GDSC.txt \
#    -out_dir results/IC50_GDSC/Normal/RGCN -nth 0 \
#    -col_cell Cell -col_drug Drug -col_ic50 LN_IC50 -cpu 4
```

### 3-2. Training Models using Whole GDSC1+2 Dataset without Splitting Data
Training models using the entire GDSC1+2 as target label dataset without splitting data is managed by the ```retrain_total.sh``` script, which sequentially executes ```train_write.sh``` and ```retrain_total.py```. The overall process is similar to that described in **Section 3-1**. Since the GDSC1+2 dataset is not split into training, validation, and test sets, early stopping is not applied when training models for 100 epochs by default.

```
bash retrain_total.sh
# python retrain_total.py \
#    -cell processed/cell_data_biocarta/SANGER_RNA_KNN5_STR9_Reg_Corr.pickle \
#    -drug processed/drug_data/GDSC_Drug_Graph.pickle \
#    -ic50 data/ic50_data/IC50_GDSC.txt \
#    -out_dir results/IC50_GDSC/Normal/RGCN \
#    -col_cell Cell -col_drug Drug -col_ic50 LN_IC50 -cpu 4 -seed_model 2021
```

## 4. Testing Models
Testing models is performed using test bash scripts (e.g., ```test_ccle.sh```, ```test_tcga.sh```, ```test_chembl.sh```), which sequentially execute ```test_write.sh``` and ```test.py```. The process is similar to **Section 3-1**. To output only predicted response values without calculating performance metrics, set the parameter ```-col_ic50``` to ```0```. We enhanced model interpretability with Grad-CAM.

### test_XXX.py
* [```-dir_param```] Model weight parameters (input, pth)
* [```-dir_hparam```] Model hyperparameters (input, Pickle)
* [```-out_file```] Prediction results (output, CSV)
* [```-out_time```] Logs for inference time when using GPU (optional output, CSV, Default: ```None```)
* [```-out_grad_cam```] Pathway importance scores with Grad-CAM (optional output, CSV, Default: ```None```)

### File Description
* ```test_ccle.sh``` : Predict IC<sub>50</sub> values in CCLE (toy example)
* ```test_tcga.sh``` : Predict drug responses in TCGA (clinical application)
* ```test_miss.sh``` : Predict missing IC<sub>50</sub> values in GDSC (case study)
* ```test_rest.sh``` : Predict IC<sub>50</sub> values for cell lines listed in SANGER Cell Passports but not estimated in GDSC (case study)
* ```test_chembl.sh``` : Predict IC<sub>50</sub> values in ChEMBL (external benchmark test)
* ```test_liu_2024_in_vivo.sh``` : Predict drug responses for Liu et al. (2024) (clinical application)

```
bash test_ccle.sh
# python test.py \
#    -cell processed/cell_data_biocarta/CCLE_RNA_KNN5_STR9_Reg_Corr.pickle \
#    -drug processed/drug_data/CCLE_Drug_Graph.pickle \
#    -ic50 data/ic50_data/IC50_CCLE.txt \
#    -dir_param results/IC50_GDSC/Normal/RGCN/param_retrain_seed2021.pt \
#    -dir_hparam results/IC50_GDSC/Normal/RGCN/hyper_param_retrain_seed2021.pickle \
#    -out_file results/IC50_GDSC/Normal/RGCN/pred_ccle_seed2021.csv \
#    -out_grad_cam results/IC50_GDSC/Normal/RGCN/gcam_ccle_seed2021.csv \
#    -col_cell Cell_BROAD_ID -col_drug Drug_CID -col_ic50 LN_IC50 -cpu 4
```
