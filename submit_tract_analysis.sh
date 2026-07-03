#!/bin/bash
#SBATCH --job-name="tract_analysis"
#SBATCH --time=0:10:00

#SBATCH --output='./logs/slurmJob_tract_analysis.out'
#SBATCH --error='./logs/slurmJob_tract_analysis.err'

json_file="INSERT/study/subjects/subj_list.json"
root="INSERT/study/"

# Read the JSON, remove brackets, and count entries
raw_json=$(<"$json_file")
num_patients=$(echo "$raw_json" | sed 's/[][]//g' | tr ',' '\n' | tr -d '"' | tr -d ' ' | wc -l)
max_index=$((num_patients - 1))

jobstr0=$(sbatch --parsable --array=0-${max_index} submit_tracking.sh "$root" "$json_file")
#jobstr0=$(sbatch --parsable --array=0-${max_index} skip.sh)
jobstr1=$(sbatch --dependency=afterany:${jobstr0##* } --array=0-${max_index} --parsable submit_unravel.sh "$root" "$json_file")
jobstr2=$(sbatch --dependency=afterany:${jobstr1##* } submit_merge_unravel.sh "$root")
sbatch --dependency=afterany:${jobstr2##* } submit_to_excel.sh "$root"