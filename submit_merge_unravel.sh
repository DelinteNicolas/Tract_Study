#!/bin/bash
#SBATCH --job-name="merge_patients"
#
#SBATCH --ntasks=2
#SBATCH --time=0:10:00
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type='FAIL'
#SBATCH --mail-user='nicolas.delinte@uclouvain.be'

#SBATCH --output='./logs/slurmJob_merge.out'
#SBATCH --error='./logs/slurmJob_merge.err'

root="$1"

python output_unravel.py "$root" "merge"