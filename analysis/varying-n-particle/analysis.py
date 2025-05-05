#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import MDAnalysis as mda
from nmrdfrommd import NMRD
from concurrent.futures import ProcessPoolExecutor, as_completed

from utilities import save_result, get_git_repo_path

def process(npart, n, git_path):
    """Process one temperature point"""
    data_dir = os.path.join(git_path, "data", "varying-n-particle", f"n{npart}")
    topology_file = os.path.join(data_dir, "topology.data")
    trajectory_file = os.path.join(data_dir, "dump"+str(n)+".xtc")

    if not os.path.exists(topology_file) or not os.path.exists(trajectory_file):
        print(f"[n={npart}] Missing MD data")
        return

    try:
        u = mda.Universe(topology_file, trajectory_file)
        all_atoms = u.select_atoms("all")

        nmr = NMRD(
            u=u,
            atom_group=all_atoms,
            number_i=1)
        nmr.run_analysis()
        
        save_result(nmr, n, name=f"n{npart}")
    except Exception as e:
        print(f"[n={npart}, n ={n}] Error: {e}")

def main(max_iterations):
    git_path = get_git_repo_path()
    n_parts = ["49", "90", "162", "292", "526", "949", "1709", "3080", "5550", "10000"]

    all_N = np.arange(1, 11)
    for iteration in range(max_iterations):
        print(f"\n--- Iteration {iteration + 1} ---")
        with ProcessPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(process, npart, n, git_path)
                       for npart in n_parts
                       for n in all_N]
            for future in as_completed(futures):
                future.result()

if __name__ == "__main__":
    main(max_iterations=3000)
