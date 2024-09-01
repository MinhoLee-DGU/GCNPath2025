#/usr/bin/env Rscript

# Even in the canonical SMILES,
# GDSC drugs are represented by ALL CAPITAL characters unlike ChEMBL
# We decided to transform canonical SIMLES from ChEMBL into those from PubChem
# Ex. [GDSC, CID 864] C1CSSC1CCCCC(=O)O
# Ex. [ChEMBL, CHEMBL137635] CN(c1ccccc1)c1ncnc2ccc(N/N=N/Cc3ccccn3)cc12

file = "_data/SMILES_ChEMBL.csv"
SMILES_ChEMBL = read.csv(file)

file = "_data/Query_ChEMBL.txt"
cat(SMILES_ChEMBL$Molecule_ChEMBL_ID, file=file, sep="\n")

file = "_data/Query_SMILES.txt"
cat(SMILES_ChEMBL$Canonical_SMILES, file=file, sep="\n")

# PubChem Identifier Exchange Service
# https://pubchem.ncbi.nlm.nih.gov/idexchange
# [1] Input ID List : SMILES
# [2] Input ID List : Registry, ChEMBL
# Output IDs : SMILES

file = "_data/Reply_ChEMBL.txt"
SMILES_ChEMBL_PCh1 = read.csv(file, sep="\t", header=F, na.strings="")
colnames(SMILES_ChEMBL_PCh1) = c("IDs_ChEMBL", "SMILES_PubChem")
SMILES_ChEMBL_PCh1$SMILES_PubChem %>% is.na %>% sum   # 8

file = "_data/Reply_SMILES.txt"
SMILES_ChEMBL_PCh2 = read.csv(file, sep="\t", header=F, na.strings="")
colnames(SMILES_ChEMBL_PCh2) = c("SMILES_ChEMBL", "SMILES_PubChem")
SMILES_ChEMBL_PCh2$SMILES_PubChem %>% is.na %>% sum   # 3876

idx1 = match(SMILES_ChEMBL$Molecule_ChEMBL_ID, SMILES_ChEMBL_PCh1$IDs_ChEMBL)
idx2 = match(SMILES_ChEMBL$Canonical_SMILES, SMILES_ChEMBL_PCh2$SMILES_ChEMBL)
smiles = coalesce(SMILES_ChEMBL_PCh1$SMILES_PubChem[idx1], SMILES_ChEMBL_PCh2$SMILES_PubChem[idx2])

SMILES_ChEMBL$Canonical_SMILES = smiles
SMILES_ChEMBL$Canonical_SMILES %>% is.na %>% sum   # 8
SMILES_ChEMBL = SMILES_ChEMBL %>% subset(!is.na(Canonical_SMILES))

file = "_data/SMILES_ChEMBL_tCNNS.csv"
write.csv(SMILES_ChEMBL, file=file, row.names=F)
