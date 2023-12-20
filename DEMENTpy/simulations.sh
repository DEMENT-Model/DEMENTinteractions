# Usage: bash simulations.sh

# Loop through each of the three different climate scenarios
for i in input/ambient/ input/drought/ input/moist
do
	# Loop through 10 different microbial initializations, using random number seeds
	for j in 12089 22765 26152 3365 4325 440 5151 5192 5732 8104
	do
		# Submit JOB-fullsims.sh to run full simulation
		export i; export j
		sbatch JOB-fullsims.sh
		
		# Loop through 1-25 taxa and conduct exclusion experiments
		for taxa in {0..24}
		do
			export taxa
		 	sbatch JOB-exclusions.sh
		done
		
	done
	
done


# Submit job to collate outputs
sbatch JOB-outputs.sh