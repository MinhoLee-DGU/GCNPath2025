#!/usr/bin/bash

data=$1
dir_out=$2
choice=$3

nth=$4
ic50_th=$5

seed=$6
retrain=$7
use_slurm=$8

CPU=20
node=full

case $ic50_th in
    0) RAM=70 ;;
    1) RAM=55 ;;
    2) RAM=40 ;;
esac

case $choice in
    0) ch=N ;;
    1) ch=C ;;
    2) ch=D ;;
    *) ch=S ;;
esac

train_file=train.py
train_file2=main.R
f_name=MF${ic50_th}_${ch}${nth}
f_name_1=exe/${f_name}.sh

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -c $CPU" >> $f_name_1
echo "#SBATCH -p $node" >> $f_name_1

echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=${RAM}GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

echo "dir_conda=\`conda info --base\`" >> $f_name_1
echo "export PATH=\$dir_conda/envs/GCNPath/bin" >> $f_name_1
echo "export PYTHONNOUSERSITE=1" >> $f_name_1

echo "python $train_file -ic50 $data -out_dir $dir_out \
    -choice $choice -nth $nth -seed_model $seed" >> $f_name_1
echo "Rscript $train_file2 $nth $dir_out $CPU" >> $f_name_1
chmod 744 $f_name_1

case $use_slurm in
    0) bash $f_name_1 ;;
    *) sbatch $f_name_1 ;;
esac
