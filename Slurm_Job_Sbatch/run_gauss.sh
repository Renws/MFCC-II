#!/bin/bash

#SBATCH --account=torch_pr_281_general

#SBATCH -J gauss
#SBATCH -o job.%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=25GB
#SBATCH -p cs
#SBATCH --time=72:00:00

source ~/.bashrc
export GAUSS_SCRDIR=$SLURM_TMPDIR
export PGI_FASTMATH_CPU=sandybridge

echo "Running on node: $(hostname)"


