#!/usr/bin/bash

ic50=$1
cell=$2
drug=$3

dir_param=$4
dir_test=$5

col_cell=$6
col_drug=$7
col_ic50=$8

RAM=8
CPU=4

f_name=$9
f_name_1=exe/$f_name.sh

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -c $CPU" >> $f_name_1
# echo "#SBATCH -p lys" >> $f_name_1
echo "#SBATCH -p glu" >> $f_name_1
# echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=${RAM}GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

echo "dir_conda=\`conda info --base\`" >> $f_name_1
echo "export PATH=\$dir_conda/envs/HiDRA/bin" >> $f_name_1

echo "python HiDRA_predict.py -ic50 $ic50 -dir_in $cell -drug $drug \
    -dir_param $dir_param -dir_test $dir_test -cpu $CPU \
    -col_cell $col_cell -col_drug $col_drug -col_ic50 $col_ic50" >> $f_name_1

chmod 777 $f_name_1
# bash $f_name_1
sbatch $f_name_1
