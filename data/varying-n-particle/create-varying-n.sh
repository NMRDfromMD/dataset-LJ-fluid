#!/bin/bash

set -e

n_part=2048
cut_off=4.0
dumping=10
T_star=0.8

for n_part in 40 61 95 148 230 356 551 854 1322 2048
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
            newline='cut_off_lj = '${cut_off}
            oldline=$(cat lj-to-real.py | grep 'cut_off_lj = ')
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
