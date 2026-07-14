#!/bin/bash

set -e

# Change accordingly
lmp=/home/simon/Softwares/lammps-22Jul2025/src/lmp_mpi

for T_star in 3.4 3.8
do
    folder=T${T_star}/
    cd $folder
        mpirun -np 8 ${lmp} -in input.lmp
    cd .. 
done
