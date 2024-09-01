#!/usr/bin/bash

nth_list=$(seq 4 30)
for nth in $nth_list
do
    bash process_assay_write.sh $nth
done
