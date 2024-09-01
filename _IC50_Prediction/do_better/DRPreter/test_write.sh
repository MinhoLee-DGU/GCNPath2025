#!/usr/bin/bash

ic50=$1
cell=$2
drug=$3

dir_param=$4
dir_test=$5

col_cell=$6
col_drug=$7
col_ic50=$8

f_name=$9
f_name_1=exe/$f_name.sh

gpu=${10}
case $gpu in
    0) CPU=2 ;;
    1) CPU=3 ;;
esac

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -n 1" >> $f_name_1
echo "#SBATCH -c $CPU" >> $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1

if [[ $gpu -eq 1 ]] ; then
	echo "#SBATCH -p glu" >> $f_name_1
    echo "#SBATCH --gres=gpu:1" >> $f_name_1
else
	echo "#SBATCH -p fMet,met" >> $f_name_1
fi

echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=4GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

echo "dir_conda=\`conda info --base\`" >> $f_name_1
echo "export PATH=\$dir_conda/envs/geometric/bin" >> $f_name_1

echo "python main.py -ic50 $ic50 -pretrain 0 --mode 'test' --sim 0 \
    -cell $cell -drug $drug -dir_param $dir_param -dir_test $dir_test \
    -col_cell $col_cell -col_drug $col_drug -col_ic50 $col_ic50 -cpu $((CPU-1))" >> $f_name_1
chmod 777 $f_name_1
sbatch $f_name_1
