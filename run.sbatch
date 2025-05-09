#!/bin/bash
#SBATCH -A m3845_g       
#SBATCH -C gpu   
#SBATCH -q regular
#SBATCH -t 20:30:00
#SBATCH -c 32          
#SBATCH --gpus-per-task=1 
#SBATCH -n 8   
#SBATCH -o output.log
#SBATCH -e error.log
module load python
srun -n 8 python main_code.py
