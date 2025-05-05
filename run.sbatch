#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 10:00:00
#SBATCH -J fast_py
#SBATCH -o output.log
#SBATCH -e error.log

module load python
srun -n 8 python main_code.py