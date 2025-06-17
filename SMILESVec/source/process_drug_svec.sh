#!/usr/bin/bash
export PYTHONNOUSERSITE=1

# Embeddings can be downloaded from the following source
# https://cmpe.boun.edu.tr/~hakime.ozturk/smilesvec.html

smi=../../data/drug_data/SMILES_GDSC.csv
out=../../processed/drug_data/GDSC_Drug_SMILESVec_C23Word.pickle
# out=../../processed/drug_data/GDSC_Drug_SMILESVec_PChWord.pickle
# out=../../processed/drug_data/GDSC_Drug_SMILESVec_Ch23Char.pickle
# out=../../processed/drug_data/GDSC_Drug_SMILESVec_PChChar.pickle

emb=../embedding/drug.l8.chembl23.canon.ws20.txt
# emb=../embedding/drug.pubchem.canon.l8.ws20.txt
# emb=../embedding/drug.chembl.canon.l1.ws20.txt
# emb=../embedding/drug.pubchem.canon.l1.ws20.txt

python getsmilesvec.py $emb $smi $out
