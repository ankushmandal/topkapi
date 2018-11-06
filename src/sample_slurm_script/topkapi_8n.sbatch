#!/bin/bash
#SBATCH --job-name=topkapi_8n
#SBATCH --output="slurm-topkapi_8n-%j.out"
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --time=00:03:00
#SBATCH --export=ALL
srun -n 8 -c 8 ./topkapi -t 8 file1 file2 file3 file4 file5 file6 file7 file8 -o my_output.txt
