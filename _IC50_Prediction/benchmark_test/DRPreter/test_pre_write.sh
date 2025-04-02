#!/usr/bin/bash

ic50=$1
seed=$2

case $ic50 in 
	0) IC50=IC50_GDSC.txt ;;
	1) IC50=IC50_GDSC1.txt ;;
	2) IC50=IC50_GDSC2.txt ;;
esac

f_name=DRP_Pre${ic50}_S${seed}
f_name_1=exe/${f_name}.sh
echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -n 1" >> $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=4GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1
echo "python main.py -ic50 $IC50 --mode 'test' -pretrain 1 --seed $seed" >> $f_name_1
chmod 777 $f_name_1
sbatch $f_name_1
