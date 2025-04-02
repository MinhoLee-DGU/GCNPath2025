# import numpy and pandas
import numpy as np
import pandas as pd
import openpyxl
import pickle
import csv
import math
import os


#Loading cell line expression
#Expression file is 'RMA normalised expression data for cell-lines' 
#from GDSC data portal(https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html)
expression_df=pd.read_csv('./_data/Cell_line_RMA_proc_basalExp.txt',sep='\t',index_col=0)
expression_df=expression_df[expression_df.columns[1:]]


#Loading cell line information
#Cell line information file is 'Annotated list of cell-lines'
#from GDSC data portal(https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html)
cellline_information=pd.read_excel('./_data/TableS1E.xlsx',index_col=0, engine='openpyxl')
cellline_information=cellline_information.iloc[3:-1]
cellline_information=cellline_information[cellline_information.columns[0:2]]
cellline_information.columns=['Cell line name','COSMIC identifier']

#Excluding cell lines whose expression values are not valid
cellline_list=cellline_information['COSMIC identifier']
cellline_list=[str(x) for x in cellline_list]
cosmic_list=expression_df.columns
cosmic_list=[x[5:] for x in cosmic_list]
isin_list=[(cosmic in cellline_list) for cosmic in cosmic_list]
expression_df=expression_df.loc[:,isin_list]   # 968

#Excluding expressions that are not gene
expression_list=expression_df.index
expression_list=[str(x) for x in expression_list]
isin_list=[x!='nan' for x in expression_list]
expression_df=expression_df.loc[isin_list]

# #Converting COSMIC identifier into Cell line name
# cellline_name_dic={}
# for idx,x in cellline_information.iterrows():
#     cellline_name_dic[str(x['COSMIC identifier'])]=x['Cell line name']
# cosmic_list=expression_df.columns
# cosmic_list=[x[5:] for x in cosmic_list]
# cellline_new_col=[cellline_name_dic[cosmic] for cosmic in cosmic_list]
# expression_df.columns=cellline_new_col

# Instead of cell names, we utilize COSMIC ID
cosmic_list=expression_df.columns
cosmic_list=[x[5:] for x in cosmic_list]
expression_df.columns = cosmic_list
isin_list = ["." not in _ for _ in expression_df.columns]
expression_df=expression_df.loc[:,isin_list]   # 968

# Save Temp RNA file...
expression_df.transpose().to_csv("_data/GDSC_RNA_Array.csv", header=True, index=True)

#Transform expression values into z-score
from scipy.stats import zscore
expression_df=expression_df.apply(zscore)

expression_df.index=expression_df.index.rename('Gene_Symbol')
# expression_df.to_csv('./expression.csv')
# Notice that z-score is performed by cells, not by genes! [Genes x Cells]
# expression_df.apply(np.mean, 0)   # Mostly zeros
# expression_df.apply(np.mean, 1)   # Not zeros


#Make gene expressions grouped into gene sets
#Loading Gene Set
#Gene Set File (gmt) is 'KEGG subset of CP' from MSigDB (http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
GeneSet_List=[]
GeneSetFile='./Training/geneset.gmt'

with open(GeneSetFile) as f:
    reader = csv.reader(f)
    data = list(list(rec) for rec in csv.reader(f, delimiter='\t')) #reads csv into a list of lists
    for row in data:
        GeneSet_List.append(row)

# GeneSet_Dic={}
pathways = list()
for GeneSet in GeneSet_List:
    pathways.append(GeneSet[0])
    # GeneSet_Dic[GeneSet[0]]=GeneSet[2:]

# Gene list are already exists...
GeneSet_Dic = {}
dir_geneset = "./Prediction/input_dir"

for i in range(186) :
    file = "{}/{}.csv".format(dir_geneset, i)
    exp_example = pd.read_csv(file, index_col=0)
    GeneSet_Dic[pathways[i]] = list(exp_example.columns.values)

geneset_int = list(set(i) for i in GeneSet_Dic.values())
print("# of Genes : {}".format(len(set.union(*geneset_int))))


#Delete genes that are not valid
#In here, E3 is just a name of one of cell line that is valid
GeneSet_Dic_withoutNA={}
key1 = set(expression_df['906826'].dropna().index)
# key1 = set(expression_df['ES3'].dropna().index)

for GeneSet in GeneSet_Dic:
    key2 = set(GeneSet_Dic[GeneSet])
    GeneSet_Dic_withoutNA[GeneSet] = list(key1 & key2)
    # GeneSet_Dic_withoutNA[GeneSet]=expression_df['ES3'][GeneSet_Dic[GeneSet]].dropna().index.values
    # KeyError: Passing list-likes to .loc or [] with any missing labels is no longer supported

expression_df=expression_df.transpose()
len(expression_df.index)        # 968 [1014]
expression_df.index.nunique()   # 968 [1014]

def CelllineFeatureExtract(ExpressionMatrix, GeneSetDic, CellLine):
    X_Feature=[]
    for GeneSet in GeneSetDic.keys():
        Gene_in_GeneSet=[]
        for Gene in GeneSetDic[GeneSet]:
            Gene_in_GeneSet.append(Gene)
        X_Feature.append(ExpressionMatrix[Gene_in_GeneSet].loc[[CellLine]])
    
    return X_Feature

cellline_input=[]
for i in range(len(GeneSet_Dic_withoutNA)):
    cellline_input.append(pd.DataFrame())
for cellline in expression_df.index:
    x=CelllineFeatureExtract(expression_df,GeneSet_Dic_withoutNA,cellline)
    for j in range(len(GeneSet_Dic_withoutNA)):
        cellline_input[j]=cellline_input[j].append(x[j])


# The Number of KEGG Genes in GDSC Array Data
# geneset_int = list(set(cellline_input[i].columns) for i in range(len(cellline_input)))
geneset_int = list(set(i) for i in GeneSet_Dic_withoutNA.values())
print("# of Genes : {}".format(len(set.union(*geneset_int))))

dir_out = "../../do_better/HiDRA/_data"
dir_exp = "{}/ProcessedFile".format(dir_out)
file_gset = "{}/geneset.pickle".format(dir_out)
os.makedirs(dir_exp, exist_ok=True)

with open(file_gset, "wb") as f :
    pickle.dump(GeneSet_Dic_withoutNA, f)

for idx, df in enumerate(cellline_input):
    df.to_csv(dir_exp+"/"+str(idx)+'.csv')
