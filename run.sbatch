#!/bin/bash
#SBATCH -A m3845_g       
#SBATCH -C gpu   
#SBATCH -q regular
#SBATCH -t 00:30:00
#SBATCH -c 32          
#SBATCH --gpus-per-task=1 
#SBATCH -n 1   
#SBATCH -o output.log
#SBATCH -e error.log
module load python
srun -n 8 python main_code.py
