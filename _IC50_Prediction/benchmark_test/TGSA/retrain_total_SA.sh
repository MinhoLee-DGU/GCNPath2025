#!/usr/bin/bash

ic50=0
choice=0
gpu=1
tgsa=1
retrain=1
nth=0

seed_list=$(seq 2022 2030)
# seed_list=(42 $(seq 2022 2030))

for seed in ${seed_list[@]} 
do
    bash train_write.sh $ic50 $choice $nth $gpu $tgsa $retrain $seed
done
