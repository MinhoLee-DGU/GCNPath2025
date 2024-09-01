import os
import json
import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client

batch_size = 1000
dir_in = "../../processed_data/ic50_data/ChEMBL"
dir_out = "../../processed_data/ic50_data/ChEMBL/Temp"
os.makedirs(dir_out, exist_ok=True)

activity_total = new_client.activity
file_batch = "{}/batch_idx.json".format(dir_out)
file_in = "{}/Assay_ChEMBL.csv".format(dir_in)

Assay_ChEMBL = pd.read_csv(file_in, low_memory=False)
assay_total = Assay_ChEMBL.Assay_ChEMBL_ID.values.tolist()
num_assay = len(Assay_ChEMBL)
print("### Batch collection started...")

if os.path.exists(file_batch) :
    with open(file_batch, "r") as f :
        batch = json.load(f)
    print("### Batch collection continued from {}...".format(batch["idx"]))
else :
    batch = {"idx":0}
    with open(file_batch, "w") as f :
        json.dump(batch, f)
    print("### Batch collection initialized...")


while True :
    idx1 = batch_size * batch["idx"]
    idx2 = batch_size * (batch["idx"]+1)
    
    if idx1 > num_assay : 
        print("# Batch start {} is larger than {}...".format(idx1, num_assay))
        print("### Batch collection completed...")
        break
    idx2 = idx2 if idx2>num_assay else idx2
    file_out = "{}/IC50_ChEMBL_Temp_{}.csv".format(dir_out, batch["idx"])
    
    with open(file_batch, "w") as f :
        json.dump(batch, f)
    
    if os.path.isfile(file_out) : 
        print("# File already exists...")
    else :
        activity_temp = activity_total.filter(assay_chembl_id__in=assay_total[idx1:idx2], type="IC50")
        activity_temp = activity_temp.filter(standard_value__isnull=False)
        activity_temp = activity_temp.filter(standard_units__isnull=False)
        activity_temp = activity_temp.filter(standard_relation__isnull=False)
        IC50_ChEMBL_Temp = pd.DataFrame(activity_temp)
        IC50_ChEMBL_Temp.to_csv(file_out, index=False)

    print("# Batch {} is completed... [n={}]".format(batch["idx"], len(IC50_ChEMBL_Temp)))
    batch["idx"] += 1
