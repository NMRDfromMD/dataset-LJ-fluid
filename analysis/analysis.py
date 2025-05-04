#!/usr/bin/env python
# coding: utf-8

import os
import MDAnalysis as mda
from nmrdfrommd import NMRD
from concurrent.futures import ProcessPoolExecutor, as_completed

from utilities import save_result, get_git_repo_path


def process_temperature(T, git_path):
    """Process one temperature point"""
    data_dir = os.path.join(git_path, "data", f"T{T}")
    topology_file = os.path.join(data_dir, "topology.data")
    trajectory_file = os.path.join(data_dir, "dump.xtc")

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
        
        save_result(nmr, name=f"T{T}")
        print(f"T={T} Success")
    except Exception as e:
        print(f"[T={T}] Error: {e}")

def main(max_iterations):
    git_path = get_git_repo_path()
    temperatures = ["0.8", "1.0", "1.2", "1.5", "1.8", "2.2", "2.6", "3.0"]

    for iteration in range(max_iterations):
        print(f"\n--- Iteration {iteration + 1} ---")
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_temperature, T, git_path) for T in temperatures]
            for future in as_completed(futures):
                future.result()  # Trigger exceptions if any

if __name__ == "__main__":
    main(max_iterations=2000)
