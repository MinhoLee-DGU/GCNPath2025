#!/usr/bin/bash

ic50=$1
choice=$2
gpu=1
retrain=1

case $choice in
    3) mode=$(seq 0 24) ;;
    *) mode=$(seq 3 9) ;;
esac

for nth in ${mode[@]}
do
    bash train_write.sh $ic50 $choice $nth $gpu $retrain
done
