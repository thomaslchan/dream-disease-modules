#!/bin/bash
#
# Usage
# -----
# $ bash launch_experiments.sh ACTION_NAME
#
# where ACTION_NAME is either 'list' or 'submit' or 'run_here'

if [[ -z $1 ]]; then
    ACTION_NAME='list'
else
    ACTION_NAME=$1
fi


## Imputation method
for impute in 'zero' 'max' 'min' 'mean' 'median' 'max_local' \
                    'min_local' 'mean_local' 'median_local'
do
    export impute=$impute
    
    ## Use this line to see where you are in the loop
    echo "impute=$impute"

    ## NOTE all env vars that have been 'export'-ed will be passed along to the .slurm file

    if [[ $ACTION_NAME == 'submit' ]]; then
        ## Use this line to submit the experiment to the batch scheduler
        sbatch < do_experiment.slurm
    
    elif [[ $ACTION_NAME == 'run_here' ]]; then
        ## Use this line to just run interactively
        bash do_experiment.slurm
    fi
done


