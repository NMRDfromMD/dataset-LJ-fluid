#!/bin/bash

set -e

# Change accordingly
lmp=/home/simon/Softwares/lammps-27Jun2024/src/lmp_mpi

for dumping in 1 2 3 4 7 10 15 22 33 49
do
    folder=d${dumping}/
    cd $folder
        mpirun -np 8 ${lmp} -in input.lmp
    cd .. 
done
