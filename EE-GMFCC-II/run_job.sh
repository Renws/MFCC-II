#!/bin/bash

#SBATCH --account=torch_pr_281_general

#SBATCH -J xyztogjf
#SBATCH -o job.%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB
#SBATCH -p cs
#SBATCH --time=72:00:00

mamba activate test

echo "Running on node: $(hostname)"

python xyztogjf_charge.py
