#!/bin/bash
 
#SBATCH --job-name=spectroscopy
#SBATCH --chdir=/work/user
#SBATCH --output=/work/%u/%x-%A-%a.log
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=5G

module load foss/2020b R/4.0.4-2
module list

input_dir="$1"

output_dir="/work/$USER/$SLURM_JOB_NAME-$SLURM_JOB_ID"
mkdir -p "$output_dir"

output="$output_dir/$SLURM_ARRAY_TASK_ID.rds"


Rscript /home/user/Spectroscopy_v_6_3.R "$input_dir" "$output" "$SLURM_ARRAY_TASK_ID"
