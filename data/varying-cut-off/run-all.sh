#!/bin/bash

set -e

# Change accordingly
lmp=/home/simon/Softwares/lammps-27Jun2024/src/lmp_mpi

for cut_off in 1 1.3 1.6 1.9 2.2 2.5 2.8 3.1 3.4 3.7 4 4.3 4.6
do
    folder=co${cut_off}/
    cd $folder
        mpirun -np 8 ${lmp} -in input.lmp
    cd .. 
done
