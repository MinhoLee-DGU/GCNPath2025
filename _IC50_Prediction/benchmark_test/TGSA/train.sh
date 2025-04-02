#!/usr/bin/bash

ic50=$1
choice=$2

case $choice in 
    3) mode=$(seq 0 24) ;;
    *) mode=$(seq 0 9) ;;
esac

gpu=1
tgsa=1
retrain=1
seed=42

for nth in $mode
do
	  bash train_write.sh $ic50 $choice $nth $gpu $tgsa $retrain $seed
done
