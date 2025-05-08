import numpy as np
import pandas as pd
import os
from math import isclose, sqrt

def generate_force_constants_and_run_phonopy(
    in_path, n, alat, a_val, id_1_lines, id_2_lines, init_no, ndisps, root=""
):
    """
    Generate FORCE_CONSTANTS file and run phonopy for B2 atom_type1 - atom_type2 structure.

    Parameters:
    - in_path: Path to simulation folder (contains input.txt and NiTi_bcc.csv).
    - n: Total number of atoms.
    - alat: Lattice parameter (used for minimum image convention).
    - a_val: Reference lattice constant (for identifying neighbors).
    - id_1_lines: List of lines for type 1 atoms.
    - id_2_lines: List of lines for type 2 atoms.
    - init_no: Initial index for displacement rows.
    - ndisps: Number of displacements.
    - root: Optional root name for phonon output naming.
    """
    ide_lat = np.genfromtxt(os.path.join(in_path, 'input.txt'))[:n, 2:]

    # Distance matrices
    ideal_distances = np.zeros((n, 3, n))
    ideal_dist_sca = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                ideal_distances[i, :, j] = ide_lat[j] - ide_lat[i]
                for d in range(3):
                    if ideal_distances[i, d, j] > (alat / 2):
                        ideal_distances[i, d, j] -= alat
                    elif ideal_distances[i, d, j] <= (-alat / 2):
                        ideal_distances[i, d, j] += alat
                ideal_dist_sca[i, j] = np.linalg.norm(ideal_distances[i, :, j])

    # Neighbor distances for 5 shells
    nn_dist = [
        (sqrt(3) * a_val) / 2,
        a_val,
        sqrt(2) * a_val,
        (sqrt(11) * a_val) / 2,
        sqrt(3) * a_val
    ]

    neighbors = []
    for i in range(n):
        i_neigh = [[i]] + [[] for _ in range(5)]
        for j in range(n):
            if i != j:
                for prox in range(5):
                    if isclose(ideal_dist_sca[i, j], nn_dist[prox], rel_tol=0.001):
                        i_neigh[prox + 1].append(j)
        neighbors.append(i_neigh)

    # Setup force constants
    fc_mat = np.zeros((n, n, 3, 3))
    base_mat = np.zeros((6, 3, 3))
    nn_matrix = 0.5 * a_val * np.array([[0., 0., 0.],
                                        [1., 1., 1.],
                                        [2., 0., 0.],
                                        [0., 2., 2.],
                                        [3., 1., 1.],
                                        [2., 2., 2.]])
    fc_arr = np.array([[0, 0, 14, 14],
                       [1, 1, 2, 2],
                       [3, 4, 14, 14],
                       [5, 6, 14, 7],
                       [8, 9, 10, 11],
                       [12, 12, 13, 13]])

    # Read force constants from CSV
    fc_per_ts_df = pd.read_csv(os.path.join(in_path, 'NiTi_bcc.csv'), skiprows=2, header=None)

        
    for m in range(init_no, init_no+ndisps):

        fc = [fc_per_ts_df.iloc[m, 0], fc_per_ts_df.iloc[m, 2], fc_per_ts_df.iloc[m, 3], fc_per_ts_df.iloc[m, 4], fc_per_ts_df.iloc[m, 5], fc_per_ts_df.iloc[m, 8], 
              fc_per_ts_df.iloc[m, 9], fc_per_ts_df.iloc[m, 10], fc_per_ts_df.iloc[m, 14], fc_per_ts_df.iloc[m, 15], fc_per_ts_df.iloc[m, 16], fc_per_ts_df.iloc[m, 17], 
              fc_per_ts_df.iloc[m, 18], fc_per_ts_df.iloc[m, 19], 0.0]

        fc_2 = [fc_per_ts_df.iloc[m, 1], fc_per_ts_df.iloc[m, 2], fc_per_ts_df.iloc[m, 3], fc_per_ts_df.iloc[m, 6], fc_per_ts_df.iloc[m, 7], fc_per_ts_df.iloc[m, 11], 
              fc_per_ts_df.iloc[m, 12], fc_per_ts_df.iloc[m, 13], fc_per_ts_df.iloc[m, 14], fc_per_ts_df.iloc[m, 15], fc_per_ts_df.iloc[m, 16], fc_per_ts_df.iloc[m, 17], 
              fc_per_ts_df.iloc[m, 20], fc_per_ts_df.iloc[m, 21], 0.0]


        for fc_set, species_range in zip([fc, fc_2], [(0, id_1_lines), (id_2_lines, n)]):
            for prox in range(6):
                base_mat[prox] = np.array([
                    [fc_set[fc_arr[prox, 0]], fc_set[fc_arr[prox, 2]], fc_set[fc_arr[prox, 2]]],
                    [fc_set[fc_arr[prox, 2]], fc_set[fc_arr[prox, 1]], fc_set[fc_arr[prox, 3]]],
                    [fc_set[fc_arr[prox, 2]], fc_set[fc_arr[prox, 3]], fc_set[fc_arr[prox, 1]]]
                ])
            for i in range(*species_range):
                for prox in range(6):
                    for j in neighbors[i][prox]:
                        fc_mat[i, j] = base_mat[prox]
                        if prox in [2, 3, 4]:
                            for r in range(1, 3):
                                if isclose(abs(ideal_distances[i, r, j]), nn_matrix[prox, 0], rel_tol=1e-4):
                                    fc_mat[i, j, [r, 0], :] = fc_mat[i, j, [0, r], :]
                                    fc_mat[i, j, :, [r, 0]] = fc_mat[i, j, :, [0, r]]
                        for r in range(3):
                            if ideal_distances[i, r, j] < 0:
                                fc_mat[i, j, r, :] *= -1
                                fc_mat[i, j, :, r] *= -1

        # Write FORCE_CONSTANTS
        out_filename = os.path.join(in_path, 'FORCE_CONSTANTS')
        with open(out_filename, 'w') as file:
            print(f"{n} {n}", file=file)
            for i in range(n):
                for j in range(n):
                    print(f"{i+1} {j+1}", file=file)
                    np.savetxt(file, np.around(fc_mat[i, j], decimals=5))

        # Run Phonopy
        os.chdir(in_path)
        os.system('phonopy -p -s band.conf')
        os.system('phonopy -p -s mesh.conf')
        os.system(f'mv band.yaml band_{m}.yaml')
        os.system(f'mv band.pdf band_{m}.pdf')
        os.system('rm phonopy.yaml')

    print("All phonon calculations completed.")