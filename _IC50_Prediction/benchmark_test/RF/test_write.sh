#!/usr/bin/bash

ic50=$1
cell=$2
drug=$3

param=$4
out_file=$5
out_time=$6

col_cell=$7
col_drug=$8
col_ic50=$9

node=${10}
jname=${11}
use_slurm=${12}

CPU=4
RAM=40

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
# echo "#SBATCH -p $node" >> $f_name_1
echo "#SBATCH -p glu" >> $f_name_1

# echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=${RAM}GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

echo "dir_conda=\`conda info --base\`" >> $f_name_1
echo "export PATH=\$dir_conda/envs/GCNPath/bin" >> $f_name_1
echo "export PYTHONNOUSERSITE=1" >> $f_name_1

data="-cell $cell -drug $drug -ic50 $ic50"
col="-col_cell $col_cell -col_drug $col_drug -col_ic50 $col_ic50"
echo "python test.py $data $col -dir_param $param \
    -out_file $out_file -out_time $out_time -cpu $CPU" >> $f_name_1
chmod 777 $f_name_1

case $use_slurm in
    0) bash $f_name_1 ;;
    *) sbatch $f_name_1 ;;
esac
