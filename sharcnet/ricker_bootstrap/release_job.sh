#!/bin/bash

mkdir -p Jobs
mkdir -p Data
# rm Skynet

# make

. par_table_generator.sh

MAX=`cat par_table.txt | wc -l`

for i in `seq 2 $MAX`; do
	sbatch single_job_sharc.sh $i
	sleep 1.0
done