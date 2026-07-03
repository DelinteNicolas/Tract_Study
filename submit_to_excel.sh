#!/bin/bash
#SBATCH --job-name="dic_to_excel"
#
#SBATCH --ntasks=2
#SBATCH --time=0:30:00
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type='FAIL'
#SBATCH --mail-user='nicolas.delinte@uclouvain.be'

#SBATCH --output='./logs/slurmJob_excel.out'
#SBATCH --error='./logs/slurmJob_excel.err'

root="$1"

python output_unravel_excel.py "$root"