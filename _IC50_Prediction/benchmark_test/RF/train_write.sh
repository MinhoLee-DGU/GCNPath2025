#!/usr/bin/bash

data=$1
dir_out=$2
choice=$3

nth=$4
ic50_th=$5

seed=$6
retrain=$7
use_slurm=$8

CPU=2
case $ic50_th in
    0) RAM=15
       CPU=5 ;;
    1) RAM=12 ;;
    2) RAM=9 ;; 
esac

node=full

case $choice in
    0) ch=N ;;
    1) ch=C ;;
    2) ch=D ;;
    *) ch=S ;;
esac

case $retrain in
    0) train_file=train.py
       f_name=RF${ic50_th}_${ch}${nth} ;;
    1) train_file=retrain_total.py
       f_name=RF_S${seed} ;;
esac

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

echo "python $train_file $data -out_dir $dir_out \
    -choice $choice -nth $nth -seed_model $seed -cpu $CPU" >> $f_name_1
chmod 744 $f_name_1

case $use_slurm in
    0) bash $f_name_1 ;;
    *) sbatch $f_name_1 ;;
esac
