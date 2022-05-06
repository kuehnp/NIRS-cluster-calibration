#!/bin/bash

input_dir=$1
if [[ -z $input_dir ]]
then
  echo "usage: bash $0 input_dir bootstraps"
  exit 1 
fi

shift
bootstraps=$1
if [[ -z $bootstraps ]]
then
  echo "usage: bash $0 input_dir bootstraps"
  exit 1 
fi

# submit 1st stage
a=$(sbatch --parsable -a 1-$bootstraps /home/user/submit_Spectroscopy.sh "$input_dir")

# submit 2nd stage waiting for 1st stage
sbatch --dependency=afterok:$a submit_merge_output.sh /work/user/spectroscopy-$a $bootstraps
