#!/bin/bash
input="snap_list.txt"

while IFS= read -r j
do
    echo 'running on snapshot '$j
    sed "s/iiii/$j/g" run_pk_matter.sh > run_pk_matter_$j.sh
    
    echo "run_pk_matter_$j.sh"
    sbatch run_pk_matter_$j.sh

done < "$input"
