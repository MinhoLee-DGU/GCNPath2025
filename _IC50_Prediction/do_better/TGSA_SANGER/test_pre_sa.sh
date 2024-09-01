#!/usr/bin/bash

ic50=$1
pretrain=$2

case $ic50 in 
	0) IC50=IC50_GDSC.txt ;;
	1) IC50=IC50_GDSC1.txt ;;
	2) IC50=IC50_GDSC2.txt ;;
esac

python main_SA.py -ic50 $IC50 --mode "test" --pretrain $pretrain
