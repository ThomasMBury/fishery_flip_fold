#!/bin/bash

# # Thomas Bury
# # PhD Candidate
# # Bauch computational epidemiology research group
# # Department of Applied Mathematics
# # Faculty of Mathematics
# # University of Waterloo

#SBATCH --mem=256MB
#SBATCH --time=0-00:01:00
#SBATCH --output=job%j.out
#SBATCH --ntasks=1

sleep 1.0
touch job$SLURM_JOB_ID.txt


