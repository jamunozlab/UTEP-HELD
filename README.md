# UTEP-HELD
Utilities to Execute Pipelines - Harmonic Ensamble Lattice Dynamics

UTEP-HELD is a Python-based automation framework designed for calculating force constants and performing phonon analysis using molecular dynamics data. This tool was developed to streamline high-throughput simulations

# Features

- Reads and organizes atomic positions and forces from MD simulations
- Computes force constants using custom BVK implementations
- Writes formatted data files for post-processing
- Generates input files for Phonopy
- Automatically runs phonon calculations

# Requirements

NumPy
Phonopy
bash (for os.system commands) 

# Usage

The code loops over all combinations of:

Defect percentages: defect_percentages
Temperature variations: n_temp_variations
Lattice parameter variations: m_lattice_variations

# To Run
Update the following in config.py and default_dict.py:

simulations_path
initial_lattice_parameter
initial_temperature
ntsteps, natoms, system_size, etc.
Phonopy config paths and file names
