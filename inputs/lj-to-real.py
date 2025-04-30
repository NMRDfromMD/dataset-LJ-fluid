#!/usr/bin/env python
# coding: utf-8

# LAMMPS Input Script
# Author: Simon Gravelle
# License: GNU General Public License v3.0
# Project: NMRfromMD – https://nmrdfrommd.github.io
# Description:Python script for simulating a Lennard-Jones fluid
# This script convert the LJ parameters from Grivet (2005, 10.1063/1.1955447,
# NMR relaxation parameters of a Lennard-Jones fluid from molecular-dynamics
# simulations), into the real units system of LAMMPS

import numpy as np
from pint import UnitRegistry
ureg = UnitRegistry()

# From Grivet 2005:
# $\rho^* = 0.84$
# $T^* \in [0.8, 3.0]$
# $N_{part} = 16384$
# cut off = 4.0 sigma

rho_star = 0.84
T_star = 3.2
sigma = 3 * ureg.angstrom # A
epsilon = 0.1 * ureg.kcal / ureg.mol # kcal/mol
mass = 1.0 * ureg.g / ureg.mol # g/mol
kB = 1.987204259e-3 * ureg.kcal / ureg.mol / ureg.K # kcal/mol/K

# Conversion parameters

tconv = np.sqrt(mass*sigma**2/epsilon)
dconv = sigma
tempconv = epsilon/kB

dt = tconv*0.0025
dt = np.round(dt.to(ureg.fs),4)
n_step = 20000 # 4*np.int32(t0/dt) # factor 4 for safety

print("Timestep = "+str(dt.to(ureg.fs)))
print("Number of step = "+str(n_step))
print("Simulation duration = "+str(np.round((dt*n_step).to(ureg.ps),3)))

rho = np.round(rho_star * mass / sigma**3, 4)

print("Density = "+str(rho))

n_part = 16384
L = np.round((n_part / (rho / mass))**(1/3),4)
print("Box size = "+str(L))

cut_off = 4.0*sigma

T = np.round(tempconv * T_star,4)
print("Temperature = "+str(T))

# Update parameters.inc

f = open("parameters.inc", "w")
f.write("# LAMMPS parameters \n")
f.write("variable mass equal " + str(mass.magnitude) + " # " + str(mass.units) + "\n")
f.write("variable sigma equal " + str(sigma.magnitude) + " # " + str(sigma.units) + "\n")
f.write("variable epsilon equal " + str(epsilon.magnitude) + " # " + str(epsilon.units) + "\n")
f.write("variable dt equal " + str(dt.magnitude) + " # " + str(dt.units) + "\n")
f.write("variable n_step equal " + str(n_step) + "\n")
f.write("variable n_part equal " + str(n_part) + "\n")
f.write("variable L equal " + str(L.magnitude)  + " # " + str(L.units) + "\n")
f.write("variable cut_off equal " + str(cut_off.magnitude)  + " # " + str(cut_off.units) + "\n")
f.write("variable T equal " + str(T.magnitude)  + " # " + str(T.units) + "\n")
f.close()

