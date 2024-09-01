#!/usr/bin/bash

ic50=$1
cell=$2
drug=$3

dir_param=$4
dir_tgsa=$5
dir_test=$6

col_cell=$7
col_drug=$8
col_ic50=$9

dir_sim=${10}
dir_cidx=${11}

gpu=1
RAM=60
f_name=${12}
f_name_1=exe/$f_name.sh

case $gpu in
    0) CPU=2 ;;
    1) CPU=3 ;;
esac

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -n 1" >> $f_name_1
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

echo "python main_SA.py -ic50 $ic50 --pretrain 0 --mode 'test' \
    -cell $cell -drug $drug -dir_param $dir_param -dir_tgsa $dir_tgsa -dir_test $dir_test \
    -dir_cidx $dir_cidx -dir_sim $dir_sim \
    -col_cell $col_cell -col_drug $col_drug -col_ic50 $col_ic50" >> $f_name_1
chmod 777 $f_name_1
sbatch $f_name_1
