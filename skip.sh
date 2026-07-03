#!/bin/bash
#
#SBATCH --job-name="skip"
#
#SBATCH --ntasks=1
#SBATCH --time=0:01:00
#SBATCH --mem-per-cpu=100
#SBATCH --mail-type='FAIL'
#SBATCH --mail-user='nicolas.delinte@uclouvain.be'

#SBATCH --output='./logs/slurmJob_tracking.out'
#SBATCH --error='./logs/slurmJob_tracking.err'

sleep 1