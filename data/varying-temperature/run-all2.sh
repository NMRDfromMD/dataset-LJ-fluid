#!/bin/bash

set -e

# Change accordingly
lmp=/home/simon/Softwares/lammps-27Jun2024/src/lmp_mpi

for T_star in 1.8 2.2 2.6 3.0
do
    folder=T${T_star}/
    cd $folder
        mpirun -np 8 ${lmp} -in input.lmp
    cd .. 
done
