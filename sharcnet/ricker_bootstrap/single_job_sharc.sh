#!/bin/bash

# # Thomas Bury
# # PhD Candidate
# # Bauch computational epidemiology research group
# # Department of Applied Mathematics
# # Faculty of Mathematics
# # University of Waterloo

#SBATCH --mem=1000MB
#SBATCH --time=0-00:10:00
#SBATCH --output=Jobs/job.%N.%j.out
#SBATCH --ntasks=1


python script_fold_gen.py `cat par_table.txt | head -n $1 | tail -n 1`
python script_flip_gen.py `cat par_table.txt | head -n $1 | tail -n 1`
