#!/bin/bash

set -e

for T in 0.8 1.0 1.2 1.5 1.8 2.2 2.6 3.0
do

    folder=T${T}/
    if [ ! -d "$folder" ];
    then
        mkdir $folder
        cd $folder
            cp ../../inputs/* .
            newline='T_star = '${T}
            oldline=$(cat lj-to-real.py | grep 'T_star = ')
            sed -i '/'"$oldline"'/c\'"$newline" lj-to-real.py
            python3 lj-to-real.py
            rm lj-to-real.py
        cd ..
    else
        echo 'Folder '${folder}'exist already, skipped'
    fi
done
