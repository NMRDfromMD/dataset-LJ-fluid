#!/bin/bash

set -e

# Change accordingly
lmp=/home/simon/Softwares/lammps-27Jun2024/src/lmp_mpi

for n_part in 49 90 162 292 526 949 1709 3080 5550 10000
do
    folder=n${n_part}/
    cd $folder
        mpirun -np 8 ${lmp} -in input.lmp
    cd .. 
done
