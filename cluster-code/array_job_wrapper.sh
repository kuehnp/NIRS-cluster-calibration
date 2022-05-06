#!/bin/bash
# this is a wrapper that starts the different stages of the script

# pass the Spectra directory introduced as comment on to the script for further use
input_dir=$1
if [[ -z $input_dir ]]
then
  echo "usage: bash $0 input_dir bootstraps"
  exit 1 
fi

# shift position, now the no of bootstraps introduced as comment gets passed on to the script
shift
bootstraps=$1
if [[ -z $bootstraps ]]
then
  echo "usage: bash $0 input_dir bootstraps"
  exit 1 
fi

# submit 1st stage, stage a starts [no of bootstraps] submit scripts and passes the input directory on to them
a=$(sbatch --parsable -a 1-$bootstraps /home/user/submit_Spectroscopy.sh "$input_dir")

# submit 2nd stage once a is finished and ok, then start submit_merge script and write to working folder
sbatch --dependency=afterok:$a submit_merge_output.sh /work/user/spectroscopy-$a $bootstraps
