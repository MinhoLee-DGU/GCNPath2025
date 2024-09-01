#!/usr/bin/bash

ic50=$1
seed_list=(2 16 33 61 79 100 220 653 1004 4001)

for seed in ${seed_list[@]}
do
    bash test_pre_write.sh $ic50 $seed
done
