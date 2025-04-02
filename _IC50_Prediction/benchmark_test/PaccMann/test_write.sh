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
f_name_1=exe/$f_name.sh

RAM=4
CPU=4

gdict=data/2128_genes.pkl
smi_lan=single_pytorch_model/smiles_language
col="-col_cell $col_cell -col_drug $col_drug -col_ic50 $col_ic50 -cpu $CPU"

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -n 1" >> $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=${RAM}GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

echo "source activate paccmann_predictor" >> $f_name_1
echo "python examples/IC50/test_paccmann.py $ic50 $cell $drug $gdict $smi_lan $dir_param $dir_test $dir_hparam $col" >> $f_name_1

chmod 777 $f_name_1
sbatch $f_name_1
