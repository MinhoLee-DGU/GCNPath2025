#!/usr/bin/bash

ic50=$1

case $ic50 in 
	0) IC50=IC50_GDSC.txt ;;
	1) IC50=IC50_GDSC1.txt ;;
	2) IC50=IC50_GDSC2.txt ;;
esac

source activate GraphDRP
python training.py -ic50 $IC50 -mode "test" -pretrain 1 

