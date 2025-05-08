#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import numpy as np
import MDAnalysis as mda
from nmrdfrommd import NMRD
from concurrent.futures import ProcessPoolExecutor, as_completed

from utilities import save_result, get_git_repo_path

def process(T, n, git_path):
    """Process one temperature point"""
    data_dir = os.path.join(git_path, "data", "varying-temperature", f"T{T}")
    topology_file = os.path.join(data_dir, "topology.data")
    trajectory_file = os.path.join(data_dir, "dump"+str(n)+".xtc")

    if not os.path.exists(topology_file) or not os.path.exists(trajectory_file):
        print(f"[T={T}] Missing MD data")
        return

    try:
        u = mda.Universe(topology_file, trajectory_file)
        all_atoms = u.select_atoms("all")

        nmr = NMRD(
            u=u,
            atom_group=all_atoms,
            number_i=1)
        nmr.run_analysis()
        
        save_result(nmr, n, name=f"T{T}")
    except Exception as e:
        print(f"[T={T}, n ={n}] Error: {e}")

def main(max_iterations):
    git_path = get_git_repo_path()
    all_T = ["0.8", "1.0", "1.2", "1.5", "1.8", "2.2", "2.6", "3.0"]
    all_N = np.arange(1, 11)
    for iteration in range(max_iterations):
        print(f"\n--- Iteration {iteration + 1} ---")
        with ProcessPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(process, T, n, git_path)
                       for T in all_T
                       for n in all_N]
            for future in as_completed(futures):
                future.result()

if __name__ == "__main__":
    main(max_iterations=5000)
