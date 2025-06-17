#!/usr/bin/bash

ic50=$1
cell=$2
drug=$3

dir_param=$4
dir_hparam=$5
dir_test=$6

col_cell=$7
col_drug=$8
col_ic50=$9

f_name=${10}
f_name_1=exe/${f_name}.sh

gdict=_data/2093_genes.pkl
smi_lan=single_pytorch_model/smiles_language
col="-col_cell $col_cell -col_drug $col_drug -col_ic50 $col_ic50"

RAM=8
gpu=${11}
case $gpu in
    0) CPU=1 ;;
    1) CPU=4 ;;
esac

echo "#!/usr/bin/bash" > $f_name_1
echo "#sbatch -c $CPU" >> $f_name_1
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

# echo "ulimit -n 2000" >> $f_name_1
# echo "export OMP_NUM_THREADS=$CPU" >> $f_name_1
echo "dir_conda=\`conda info --base\`" >> $f_name_1
echo "export PATH=\$dir_conda/envs/paccmann_predictor/bin" >> $f_name_1

echo "python examples/IC50/test_paccmann.py $ic50 $cell $drug \
    $gdict $smi_lan $dir_param $dir_test $dir_hparam $col -cpu $CPU" >> $f_name_1

chmod 777 $f_name_1
sbatch $f_name_1
