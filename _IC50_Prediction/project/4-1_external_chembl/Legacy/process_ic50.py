import os
import sys
import argparse
import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client

parser = argparse.ArgumentParser()
parser.add_argument("-total", type=int, default=0)
args = parser.parse_args()

# DATA_UPLOAD_MAX_MEMORY_SIZE = 100000000
dir_ = "../../processed_data/ic50_data/ChEMBL"

if args.total==0 :
    file = "{}/Assay_ChEMBL_Sub.csv".format(dir_)
    file_out = "{}/IC50_ChEMBL_Sub_Temp.csv".format(dir_)
else :
    file = "{}/Assay_ChEMBL_Total.csv".format(dir_)
    file_out = "{}/IC50_ChEMBL_Total_Temp.csv".format(dir_)

activity_total = new_client.activity
print("File in : {}".format(file))
print("File out : {}".format(file_out))

Anno_Assay_ChEMBL = pd.read_csv(file, low_memory=False)
assay_total = Anno_Assay_ChEMBL.Assay_ChEMBL_ID.values.tolist()
activities = activity_total.filter(assay_chembl_id__in=assay_total, type="IC50")
                                   
IC50_ChEMBL = pd.DataFrame(activities)
IC50_ChEMBL.to_csv(file_out, index=False)
print("# Processing is completed... [n={}]".format(len(IC50_ChEMBL)))
