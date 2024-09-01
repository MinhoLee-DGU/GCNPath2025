#!/usr/bin/bash

ic50=0
choice=0
gpu=0
retrain=1
nth=0

seed_list=$(seq 2021 2030)

for seed in ${seed_list[@]} 
do
    bash train_write.sh $ic50 $choice $nth $gpu $retrain $seed
done
