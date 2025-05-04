#!/bin/bash

set -e

n_part=16384
cut_off=4.0
dumping=10
T_star=3.0

for n_part in 200 326 532 868 1417 2312 3772 6154 10042 16383
do
    folder=n${n_part}/
    if [ ! -d "$folder" ];
    then
        mkdir $folder
        cd $folder
            cp ../../../inputs/* .

            # Replace temperature in file
            newline='T_star = '${T_star}
            oldline=$(cat lj-to-real.py | grep 'T_star = ')
            sed -i '/'"$oldline"'/c\'"$newline" lj-to-real.py

            # Replace N particle in file
            newline='n_part = '${n_part}
            oldline=$(cat lj-to-real.py | grep 'n_part = ')
            sed -i '/'"$oldline"'/c\'"$newline" lj-to-real.py

            # Replace cut off in file
            newline='cut_off = '${n_part}'*sigma'
            oldline=$(cat lj-to-real.py | grep 'cut_off = ')
            sed -i '/'"$oldline"'/c\'"$newline" lj-to-real.py

            # Replace dumping in file
            newline='dumping = '${dumping}
            oldline=$(cat lj-to-real.py | grep 'dumping = ')
            sed -i '/'"$oldline"'/c\'"$newline" lj-to-real.py

            python3 lj-to-real.py
            rm lj-to-real.py
        cd ..
    else
        echo 'Folder '${folder}'exist already, skipped'
    fi
done
