#!/usr/bin/env python

import pickle

with open("mut_dict", "rb") as f : 
    mut_dict = pickle.load(f)

feat_names = list(mut_dict.keys())
feat_index = list(mut_dict.values())
Feat_Index = pd.DataFrame({"Feature":feat_names, "Index":feat_index})
Feat_Index.to_csv("GraphDRP_Feat.txt", index=False, sep="\t")
