from wordextract import *
from cmethods import *
import sys, pickle

import torch
import numpy as np
import pandas as pd
# smiles_path = "utils/smiles_sample.txt"
emb_path = sys.argv[1]
in_path = sys.argv[2]
out_path = sys.argv[3]

def loadEmbeddings(LRNPATH):
    embeddings_index = {}

    f = open(os.path.join(LRNPATH)) #'word.11l.100d.txt'
    next(f)
    vsize = 0
    for line in f:
        values = line.split()
        word = values[0]
        vsize = len(values)-1
        coefs = np.asarray(values[1:], dtype='float32')
        embeddings_index[word] = coefs
    f.close()

    return embeddings_index, vsize


def getSMIVector(LINGOembds, smiles, q=8, wordOrChar="wd"):


    lingoList = []
    if wordOrChar == "wd":
        lingoList = createLINGOs(smiles, q)
    #elif wordOrChar == "ch":
    #    lingoList = createCHRs(smiles, "l") #ligand, q=1

    smilesVec = vectorAddAvg(LINGOembds, lingoList)

    return smilesVec


# def returnSMIVector(emb_file, smiles_path):
def returnSMIVector(emb_file):
    EMB, vsize = loadEmbeddings(emb_file)
    drug_data = pd.read_csv(in_path, sep=",")
    drugs = drug_data["Drug_CID"]
    smiles = drug_data["SMILES_CAN"]
    # smiless = [line.strip() for line in open(smiles_path)]
    
    print("Constructing SMILES vectors..")
    # smiVectors = [] 
    drug_dict = {}
    # for smi in smiless:
    for drug, smi in zip(drugs, smiles):
        v = getSMIVector(EMB, smi)
        drug_dict[drug] = torch.tensor(v)
        # smiVectors.append(getSMIVector(EMB, smi))
    
    print(f"Dim of drugs : {drug_dict[drug].shape}")
    with open(out_path, 'wb') as (f):
        pickle.dump(drug_dict, f)
    
    # pickle.dump(smiVectors, open("smiles.vec",'wb')) 
    print("Done.")


if __name__=="__main__":
    #emb_file, smiles_path
    returnSMIVector(emb_path)
    # returnSMIVector(emb_path, smiles_path)
