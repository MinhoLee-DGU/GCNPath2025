#!/usr/bin/bash

nth=$1
f_name=IC50_ChEMBL_$nth
f_name_1=exe/$f_name.sh

echo "#!/usr/bin/bash" > $f_name_1
echo "#SBATCH -J $f_name" >> $f_name_1
echo "#SBATCH -c 1" >> $f_name_1
echo "#SBATCH -p met" >> $f_name_1

echo "#SBATCH --mem=1GB" >> $f_name_1
echo "#SBATCH -o out/%j.out" >> $f_name_1
echo "#SBATCH -e out/e%j.out" >> $f_name_1
echo "source activate Process" >> $f_name_1

echo "python process_ic50_batch.py -nth $nth -total 1" >> $f_name_1
chmod 777 $f_name_1
sbatch $f_name_1
