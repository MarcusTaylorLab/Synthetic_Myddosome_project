#!/bin/bash -l

# get the variables defined in user_parameters.sh
source user_parameters.sh

# submit the slurm script with the number of arrays corresponding to number of structures to be analysed and the slurm.log directory specified
sbatch -a 0-$NUM -D $MYROSETTA/slurm_logs relax_ddg.sh

echo "Submitted jobs"