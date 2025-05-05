#!/bin/bash

set -e

# Change accordingly
lmp=/home/simon/Softwares/lammps-27Jun2024/src/lmp_mpi

for n_part in 40 61 95 148 230 356 551 854 1322 2048
do
    folder=n${n_part}/
    cd $folder
        mpirun -np 8 ${lmp} -in input.lmp
    cd .. 
done
