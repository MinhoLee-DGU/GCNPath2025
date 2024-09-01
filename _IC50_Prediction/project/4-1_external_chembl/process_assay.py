import os
import sys
import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client

cell_total = new_client.cell_line
cell_human = cell_total.filter(cell_source_organism="Homo sapiens")   # 1483
Anno_Cells_ChEMBL = pd.DataFrame(cell_human)

dir_ = "../../processed_data/ic50_data/ChEMBL"
file = "{}/Anno_Cells_ChEMBL_Temp.csv".format(dir_)
Anno_Cells_ChEMBL.to_csv(file, index=False)


assay_total = new_client.assay
cells = Anno_Cells_ChEMBL.cell_chembl_id.values.tolist()
assays = assay_total.filter(cell_chembl_id__in=cells)
Anno_Assay_ChEMBL = pd.DataFrame(assays)

dir_ = "../../processed_data/ic50_data/ChEMBL"
file = "{}/Anno_Assay_ChEMBL.csv".format(dir_)
Anno_Assay_ChEMBL.to_csv(file, index=False)

print("# Processing is completed... [n={}]".format(len(Anno_Assay_ChEMBL)))
