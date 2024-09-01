#!/usr/bin/bash

ic50=$1
choice=$2
nth=$3

tname=$4
cell=$5
drug=$6

col_cell=$7
col_drug=$8
col_ic50=$9

gpu=${10}
seed=${11}
RAM=80

case $gpu in
    0) CPU=10 ;;
    *) CPU=4 ;;
esac

case $choice in
    0) ch="N" ;;
    1) ch="C" ;;
    2) ch="D" ;;
    *) ch="S" ;;
esac

f_name="GraphDRP_T_${ch}${nth}"
f_name_1=exe/${f_name}.sh

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -c $CPU" >> $f_name_1

if [[ $gpu -eq 1 ]] ; then
    echo "#SBATCH -p glu" >> $f_name_1
	echo "#SBATCH --gres=gpu:1" >> $f_name_1
else
	echo "#SBATCH -p leu" >> $f_name_1
fi

echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=${RAM}GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

echo "dir_conda=\`conda info --base\`" >> $f_name_1
# echo "export PATH=\$dir_conda/envs/GraphDRP/bin" >> $f_name_1
echo "export PATH=\$dir_conda/envs/geometric/bin" >> $f_name_1
echo "export PYTHONNOUSERSITE=1" >> $f_name_1

echo "python retrain_total.py -ic50 $ic50 -choice $choice -nth $nth \
    -mode "test" -pretrain 0 -test_name $tname -dir_cell $cell -dir_drug $drug \
    -col_cell $col_cell -col_drug $col_drug -col_ic50 $col_ic50 -cpu $(($CPU-1)) -seed $seed" >> $f_name_1
chmod 777 $f_name_1
# bash $f_name_1
sbatch $f_name_1
