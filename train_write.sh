#!/usr/bin/bash

data=$1
dir_out=$2
choice=$3

nth=$4
option="$5"
option_attn="$6"

gpu=$7
node=$8
ic50_th=$9
seed=${10}
retrain=${11}
use_slurm=${12}
RAM=4

case $gpu in
    0) CPU=2 ;;
    *) CPU=4 ;;
esac

case $choice in
    0) ch=N ;;
    1) ch=C ;;
    2) ch=D ;;
    *) ch=S ;;
esac

case $node in
	0) node=full ;;
	1) node=full_wt_gpu ;;
	2) node=cath ;;
	*) node=$node ;;
esac

case $retrain in
    0) train_file=train.py
       f_name=GCN${ic50_th}_${ch}${nth} ;;
    # 1) train_file=retrain.py
    #    f_name=GCN${ic50_th}_${ch}${nth}R ;;
    1) train_file=retrain_total.py
       f_name=GCN_Seed${seed} ;;
esac

dir_out_=$(echo $dir_out | rev | cut -d "/" -f1 | rev)
f_name_1=exe/${f_name}_${dir_out_}.sh

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -c $CPU" >> $f_name_1

if [[ $gpu -eq 1 ]] ; then
	echo "#SBATCH -p glu" >> $f_name_1
	echo "#SBATCH --gres=gpu:1" >> $f_name_1
else 
    echo "#SBATCH -p $node" >> $f_name_1
fi

echo "#SBATCH -x lysine" >> $f_name_1
echo "#SBATCH --mem=${RAM}GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1

echo "dir_conda=\`conda info --base\`" >> $f_name_1
echo "export PATH=\$dir_conda/envs/GCNPath/bin" >> $f_name_1
echo "export PYTHONNOUSERSITE=1" >> $f_name_1

echo "python $train_file $data $option $option_attn \
    -out_dir $dir_out -choice $choice -nth $nth \
    -cpu $(($CPU-1)) -seed_model $seed" >> $f_name_1
chmod 744 $f_name_1

case $use_slurm in
    0) bash $f_name_1 ;;
    *) sbatch $f_name_1 ;;
esac
