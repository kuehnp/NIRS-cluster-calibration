#!/bin/bash
 
#SBATCH --job-name=spectroscopy_merger
#SBATCH --chdir=/work/user
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=5G

module load foss/2020b R/4.0.4-2
module list

input_dir="$1"
res_no="$2"

Rscript /home/user/merge_output.R "$input_dir" "$res_no"