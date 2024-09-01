#!/usr/bin/bash

nth_list=$(seq 1 4)
for nth in $nth_list
do
    bash process_ic50_batch_write.sh $nth
done
