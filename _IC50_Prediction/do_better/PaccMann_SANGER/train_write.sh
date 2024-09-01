#!/usr/bin/bash

ic50=$1
choice=$2
nth=$3
gpu=$4
retrain=$5
seed=$6
RAM=60

case $gpu in
    0) CPU=2 ;;
    1) CPU=4 ;;
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
    0) file=train_paccmann.py
       f_name=Pacc_G$1_${ch}$3 ;;
    # *) file=train_paccmann_retrain.py ;;
    *) file=retrain_total.py
       f_name=Pacc_Seed${seed} ;;
esac

f_name_1=exe/${f_name}.sh

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -n 1" >> $f_name_1
echo "#SBATCH -c $CPU" >> $f_name_1

if [ $gpu -eq 1 ] ; then
	echo "#SBATCH -p glu" >> $f_name_1 
	echo "#SBATCH --gres=gpu:1" >> $f_name_1
else 
	echo "#SBATCH -p full" >> $f_name_1
fi

echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=${RAM}GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

# echo "ulimit -n 2000" >> $f_name_1
# echo "export OMP_NUM_THREADS=$CPU" >> $f_name_1
echo "dir_conda=\`conda info --base\`" >> $f_name_1
echo "export PATH=\$dir_conda/envs/paccmann_predictor/bin" >> $f_name_1

case $retrain in 
    0) echo "python examples/IC50/$file -ic50 _data/$IC50 \
          -choice $choice -nth $nth -cpu $(($CPU-1))" >> $f_name_1 ;;
    *) echo "python examples/IC50/$file -ic50 _data/$IC50 \
          -choice $choice -nth $nth -cpu $(($CPU-1)) -seed $seed" >> $f_name_1 ;;
esac

chmod 744 $f_name_1
# bash $f_name_1
sbatch $f_name_1
