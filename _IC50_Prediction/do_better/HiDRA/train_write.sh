#!/usr/bin/bash

ic50=$1
choice=$2
nth=$3
gpu=$4
retrain=$5
seed=$6
RAM=6

case $gpu in
    0) CPU=10 ;;
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

case $retrain in 
    0) file=HiDRA_training.py
       f_name=HiDRA_G$1_${ch}$3 ;;
    # *) file=HiDRA_training_retrain.py
    #    f_name=HiDRA_G$1R_${ch}$3 ;;
    *) file=retrain_total.py
       f_name=HiDRA_Seed${seed}
esac

f_name_1=exe/${f_name}.sh

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -c $CPU" >> $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1

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

echo "source activate HiDRA" >> $f_name_1
case $retrain in
    0) echo "python $file -ic50 $IC50 -choice $choice \
        -dir_in _data/ProcessedFile -nth $nth -cpu $(($CPU-1))" >> $f_name_1 ;;
    1) echo "python $file -ic50 $IC50 -choice $choice \
        -dir_in _data/ProcessedFile -nth $nth -cpu $(($CPU-1)) -seed $seed" >> $f_name_1 ;;
esac

chmod 777 $f_name_1
# bash $f_name_1
sbatch $f_name_1

