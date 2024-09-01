import os
import sys
import argparse
import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client

parser = argparse.ArgumentParser()
parser.add_argument("-nth", type=int, default=0)
parser.add_argument("-total", type=int, default=0)
args = parser.parse_args()

batch_size = 50000
# DATA_UPLOAD_MAX_MEMORY_SIZE = 100000000
dir_ = "../../processed_data/ic50_data/ChEMBL"

if args.total==0 :
    file = "{}/Assay_ChEMBL_Sub.csv".format(dir_)
    file_out = "{}/IC50_ChEMBL_Sub_Temp_{}.csv".format(dir_, args.nth)
else :
    file = "{}/Assay_ChEMBL_Total.csv".format(dir_)
    file_out = "{}/IC50_ChEMBL_Total_Temp_{}.csv".format(dir_, args.nth)

activity_total = new_client.activity
print("File in : {}".format(file))
print("File out : {}".format(file_out))

Anno_Assay_ChEMBL = pd.read_csv(file, low_memory=False)
assay_total = Anno_Assay_ChEMBL.Assay_ChEMBL_ID.values.tolist()

idx1 = batch_size * args.nth
idx2 = batch_size * (args.nth+1)
num_assay = len(Anno_Assay_ChEMBL)
print("### Batch {}th [{}, {}]".format(args.nth, idx1, idx2))

if idx1 > num_assay :
    print("# Batch start is {}, larger than {}...".format(idx1, num_assay))
    sys.exit()
idx2 = idx2 if idx2>num_assay else idx2

if os.path.isfile(file_out) :
    print("# File already exists...")
else :
    activity_temp = activity_total.filter(assay_chembl_id__in=assay_total[idx1:idx2], type="IC50")
    IC50_ChEMBL_Temp = pd.DataFrame(activity_temp)
    IC50_ChEMBL_Temp.to_csv(file_out, index=False)
    print("# Processing batch is completed... [n={}]".format(len(IC50_ChEMBL_Temp)))
