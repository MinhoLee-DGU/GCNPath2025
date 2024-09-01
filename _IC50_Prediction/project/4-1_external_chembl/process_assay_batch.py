import os
import sys
import argparse
import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client

batch_size = 10
assay_total = new_client.assay

parser = argparse.ArgumentParser()
parser.add_argument("-nth", type=int, default=0)
args = parser.parse_args()

dir_ = "../../processed_data/ic50_data/ChEMBL"
file = "{}/Anno_Cells_ChEMBL.csv".format(dir_)
Anno_Cells_ChEMBL = pd.read_csv(file)
cells = Anno_Cells_ChEMBL.cell_chembl_id.values.tolist()

idx1 = batch_size * args.nth
idx2 = batch_size * (args.nth+1)
num_cells = len(Anno_Cells_ChEMBL)
print("### Batch {}th [{}, {}]".format(args.nth, idx1, idx2))

if idx1 > num_cells :
    print("# Batch start is {}, larger than {}...".format(idx1, num_cells))
    sys.exit()
idx2 = idx2 if idx2>num_cells else idx2


dir_ = "../../processed_data/ic50_data/ChEMBL/Temp"
file = "{}/Anno_Assay_ChEMBL_{}.csv".format(dir_, args.nth)
os.makedirs(dir_, exist_ok=True)

if os.path.isfile(file) :
    print("# File already exists...")
else :
    assays_temp = assay_total.filter(cell_chembl_id__in=cells[idx1:idx2])
    Anno_Assay_Temp = pd.DataFrame(assays_temp)
    Anno_Assay_Temp.to_csv(file, index=False)
    print("# Processing batch is completed... [n={}]".format(len(Anno_Assay_Temp)))
