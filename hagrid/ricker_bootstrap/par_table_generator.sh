#!/bin/bash

rm par_table.txt
touch par_table.txt

declare -a TMAX_VALS=(500);
declare -a SEED_VALS=($(seq 0 1 100));
declare -a SIGMA_VALS=(0.04);
declare -a SPAN_VALS=(0.5);
declare -a RW_VALS=(0.4);
declare -a HAM_LENGTH_VALS=(40);
declare -a HAM_OFFSET_VALS=(0.5);
declare -a SWEEP_VALS=('true');
declare -a BLOCK_SIZE_VALS=(20);
declare -a BS_TYPE_VALS=('Stationary');
declare -a N_SAMPLES_VALS=(100);



echo "tmax seed sigma span rw ham_length ham_offset sweep block_size bs_type n_samples" >> par_table.txt;

for tmax in "${TMAX_VALS[@]}"; do
for sigma in "${SIGMA_VALS[@]}"; do
for seed in "${SEED_VALS[@]}"; do
	for span in "${SPAN_VALS[@]}"; do
	for rw in "${RW_VALS[@]}"; do
	for ham_length in "${HAM_LENGTH_VALS[@]}"; do
		for ham_offset in "${HAM_OFFSET_VALS[@]}"; do
		for sweep in "${SWEEP_VALS[@]}"; do
		for block_size in "${BLOCK_SIZE_VALS[@]}"; do
			for bs_type in "${BS_TYPE_VALS[@]}"; do
			for n_samples in "${N_SAMPLES_VALS[@]}"; do


				echo "$tmax $seed $sigma $span $rw $ham_length $ham_offset $sweep $block_size $bs_type $n_samples" >> par_table.txt;

			done
			done
		done
		done
		done
	done
	done
	done
done
done
done

