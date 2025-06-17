#!/usr/bin/bash

ic50_th=0
choice=0

seed=2021
retrain=1

# cell=_data/SANGER_RNA.pickle
cell=_data/SANGER_RNA_Noise.pickle
drug=_data/GDSC_Drug_Morgan.pickle

case $ic50_th in
    0) IC50=IC50_GDSC.txt ;;
    1) IC50=IC50_GDSC1.txt ;;
    2) IC50=IC50_GDSC2.txt ;;
esac

ic50_=(IC50_GDSC IC50_GDSC1 IC50_GDSC2)
choice_=(Normal Cell_Blind Drug_Blind Strict_Blind)
# dir_out=Results/${ic50_[$ic50_th]}/${choice_[$choice]}
dir_out=Results_Noise/${ic50_[$ic50_th]}/${choice_[$choice]}

ic50=_data/$IC50
data="-cell $cell -drug $drug -ic50 $ic50"

nth=0
use_slurm=1
seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]}
do
    bash train_write.sh "$data" $dir_out \
        $choice $nth $ic50_th $seed $retrain $use_slurm
done
