#!/bin/bash

set -e

n_part=2048
cut_off=4.0
dumping=10
T_star=0.8

for T_star in 0.8 1.0 1.2 1.5 1.8 2.2 2.6 3.0
do
    folder=T${T_star}/
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
