#!/usr/bin/bash

ic50=$1
choice=$2
nth=$3
gpu=$4
tgsa=$5
retrain=$6
seed=$7
RAM=60

case $gpu in
    0) CPU=1 ;;
    *) CPU=4 ;;
esac

case $ic50 in 
	0) IC50="IC50_GDSC.txt" ;;
	1) IC50="IC50_GDSC1.txt" ;;
	2) IC50="IC50_GDSC2.txt" ;;
esac

case $choice in
	0) ch="N" ;;
	1) ch="C" ;;
    2) ch="D" ;;
    *) ch="S" ;;
esac

case $tgsa+$retrain in 
    0+0) file=main.py 
         f_name=TGDRP_G$1_${ch}$3 ;;
    # 0+1) file=main_retrain.py
    #      f_name=TGDRP_G$1R_${ch}$3 ;;
    0+1) file=retrain_total.py
         f_name=TGDRP_Seed${seed} ;;
    1+0) file=main_SA.py
         f_name=TGSA_G$1_${ch}$3 ;;
    # 1+1) file=main_SA_retrain.py
    #      f_name=TGSA_G$1R_${ch}$3 ;;
    1+1) file=retrain_total_SA.py
         f_name=TGSA_Seed${seed} ;;
esac

f_name_1=exe/${f_name}.sh

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -c $CPU" >> $f_name_1

if [[ $gpu -eq 1 ]] ; then
	echo "#SBATCH -p glu" >> $f_name_1
    echo "#SBATCH --gres=gpu:1" >> $f_name_1
else
	echo "#SBATCH -p full" >> $f_name_1
fi

echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=${RAM}GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

echo "dir_conda=\`conda info --base\`" >> $f_name_1
echo "export PATH=\$dir_conda/envs/geometric/bin" >> $f_name_1
echo "python $file -ic50 $IC50 -choice $choice \
    -nth $nth --mode train -cpu $(($CPU-1)) --seed $seed" >> $f_name_1
chmod 777 $f_name_1
# bash $f_name_1
sbatch $f_name_1
