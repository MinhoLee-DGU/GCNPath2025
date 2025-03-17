# **File Structure – ```project```**
* ```1-1_process_rna``` : Processes cell TPM data from SANGER Cell Model Passports
* ```1-2_process_net``` : Processes gene-level networks from STRING and RegNetwork
* ```1-3_gsva_ppa``` : Creates pathway graphs and examines batch correction via GSVA
* ```1-4_process_mut_cnv``` : Processes mutation and CNV data from SANGER Cell Model Passports for training TGSA
* ```2-1_get_drug_cid``` : Retrieves PubChem CIDs and SMILES for drugs in GDSC
* ```3-1_process_ic50``` : Processes IC50 data from GDSC, averaging duplicates
* ```3-2_process_ic50_ccle``` : Processes IC50 data from CCLE
* ```4-1_external_chembl``` : Processes drug and IC50 data from ChEMBL
* ```4-2_tcga_response``` : Retrieves clinical and TPM data from TCGA
* ```prepare[Model]``` : Processes cell omics data from GDSC or CCLE DepMap for training and testing each model
* ```prepare[Model]_sanger``` : Processes cell omics data from SANGER Cell Model Passports for training and testing each model
* ```functions.R``` : List of functions for processing data in R
* ```utils_gsva.R``` : Functions for compressing cell data to the pathway level using GSVA
