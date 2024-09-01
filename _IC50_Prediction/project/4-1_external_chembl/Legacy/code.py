import os
import sys
import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client

prepared_cell = True
dir_ = "../../processed_data/ic50_data/ChEMBL"
file = "{}/Anno_Cells_ChEMBL.csv".format(dir_)

if not prepared_cell :
    cell_total = new_client.cell_line
    assay_total = new_client.assay
    activity_total = new_client.activity
    cell_human = cell_total.filter(cell_source_organism="Homo sapiens")   # 1483
    
    Anno_Cells_ChEMBL = pd.DataFrame(cell_human)
    Anno_Cells_ChEMBL.to_csv(file, index=False)
else :
    Anno_Cells_ChEMBL = pd.read_csv(file)


prepared_assay = False
dir_ = "../../processed_data/ic50_data/ChEMBL/Temp"
files = os.listdir(path=dir_)
files = [os.path.join(dir_, f) for f in files]

if not prepared_assay :
    print("Process assay information first...")
    sys.exit()
else :
    for i in range(len(files)) :
        Temp = pd.read_csv(files[i])
        if i==0 :
            Anno_Assay_ChEMBL = Temp
        else : 
            Anno_Assay_ChEMBL = pd.concat([Anno_Assay_ChEMBL, Temp])


# ### Example
# # Cell      CHEMBL3308068
# # Assay     CHEMBL1648832
# ex_cell = cell_total.filter(cell_chembl_id="CHEMBL3308068")            # 1
# ex_assay = assay_total.filter(cell_chembl_id="CHEMBL3308068")          # 20
# ex_activity = activity_total.filter(assay_chembl_id="CHEMBL1648832")   # 2
