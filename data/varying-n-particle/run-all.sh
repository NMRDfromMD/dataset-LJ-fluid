#!/bin/bash

set -e

# Change accordingly
lmp=/home/simon/Softwares/lammps-27Jun2024/src/lmp_mpi

for n_part in 200 326 532 868 1417 2312 3772 6154 10042 16383
do
    folder=n${n_part}/
    cd $folder
        mpirun -np 8 ${lmp} -in input.lmp
    cd .. 
done
