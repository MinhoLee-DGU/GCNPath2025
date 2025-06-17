#!/usr/bin/bash

ic50=$1
cell=$2
drug=$3

param=$4
hparam=$5
out_file=$6
out_cam="$7"

col_cell=$8
col_drug=$9
col_ic50=${10}

gpu=${11}
node=${12}
jname=${13}
use_slurm=${14}
RAM=4

case $gpu in
    0) CPU=1 ;;
    *) CPU=4 ;;
esac

case $node in
	0) node="full" ;;
	1) node="full_wt_gpu" ;;
	2) node="cath" ;;
	*) node=$node ;;
esac

f_name=$jname
f_name_1=exe/${jname}.sh

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -c $CPU" >> $f_name_1

if [[ $gpu -eq 1 ]] ; then
	echo "#SBATCH -p glu" >> $f_name_1
	echo "#SBATCH --gres=gpu:1" >> $f_name_1
else 
    echo "#SBATCH -p $node" >> $f_name_1
fi

# echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=${RAM}GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

echo "dir_conda=\`conda info --base\`" >> $f_name_1
echo "export PATH=\$dir_conda/envs/GCNPath/bin" >> $f_name_1
echo "export PYTHONNOUSERSITE=1" >> $f_name_1

if [[ $gpu -eq 0 ]] ; then
    echo "export CUDA_VISIBLE_DEVICES=none" >> $f_name_1
fi

data="-cell $cell -drug $drug -ic50 $ic50"
col="-col_cell $col_cell -col_drug $col_drug -col_ic50 $col_ic50"
echo "python test.py $data $col -dir_param $param -dir_hparam $hparam \
    -out_file $out_file -out_grad_cam $out_cam -cpu $CPU" >> $f_name_1
chmod 777 $f_name_1

case $use_slurm in
    0) bash $f_name_1 ;;
    *) sbatch $f_name_1 ;;
esac

