#!/bin/bash

#SBATCH --account=torch_pr_281_general

#SBATCH -J 5vfk
#SBATCH -o job.%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB
#SBATCH -p cs
#SBATCH --time=72:00:00

source ~/.bashrc
mamba activate test

echo "Running on node: $(hostname)"

for i in `cat raw`
do
	./sum_ene.sh 5vfk_${i}
done > 5vfk_time_results.txt

