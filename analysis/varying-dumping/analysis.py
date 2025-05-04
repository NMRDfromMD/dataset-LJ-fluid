#!/usr/bin/env python
# coding: utf-8

import os
import MDAnalysis as mda
from nmrdfrommd import NMRD
from concurrent.futures import ProcessPoolExecutor, as_completed

from utilities import save_result, get_git_repo_path

def process(n, git_path):
    """Process one temperature point"""
    data_dir = os.path.join(git_path, "data", "varying-n-particle", f"n{n}")
    topology_file = os.path.join(data_dir, "topology.data")
    trajectory_file = os.path.join(data_dir, "dump.xtc")

    if not os.path.exists(topology_file) or not os.path.exists(trajectory_file):
        print(f"[n={n}] Missing MD data")
        return

    try:
        u = mda.Universe(topology_file, trajectory_file)
        all_atoms = u.select_atoms("all")

        nmr = NMRD(
            u=u,
            atom_group=all_atoms,
            number_i=1)
        nmr.run_analysis()
        
        save_result(nmr, name=f"n{n}")
        print(f"n={n} Success")
    except Exception as e:
        print(f"[n={n}] Error: {e}")

def main(max_iterations):
    git_path = get_git_repo_path()
    n_parts = ["200", "326", "532", "868", "1417", "2312", "3772", "6154", "10042", "16383"]

    for iteration in range(max_iterations):
        print(f"\n--- Iteration {iteration + 1} ---")
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process, n, git_path) for n in n_parts]
            for future in as_completed(futures):
                future.result()

if __name__ == "__main__":
    main(max_iterations=3000)
