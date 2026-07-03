#!/bin/bash
#SBATCH --job-name="unravel"
#
#SBATCH --ntasks=2
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type='FAIL'

#SBATCH --output='./logs/slurmJob_unravel.out'
#SBATCH --error='./logs/slurmJob_unravel.err'

root="$1"
json_file="$2"

# Process JSON to extract patient list
raw_json=$(<"$json_file")
cleaned_json=$(echo "$raw_json" | sed 's/[][]//g' | tr ',' '\n' | tr -d '"' | tr -d ' ')
mapfile -t patients <<< "$cleaned_json"

python output_unravel.py "$root" ${patients[${SLURM_ARRAY_TASK_ID}]}