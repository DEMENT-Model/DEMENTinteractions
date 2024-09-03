# Usage: bash simulations.sh

# Loop through each of the three different climate scenarios
for i in input/ambient/ input/drought/ input/moist
do
	# Loop through 10 different microbial initializations, using random number seeds
	for j in 12089 22765 26152 3365 4325 440 5151 5192 5732 8104 13089 22775 25152 335 4305 410 5051 5792 572 8141 1089 2275 2512 35 4105 100 5451 9200 5882 8417 11089 27765 29152 4365 432 7740 1510 1923 5932 1014 1209 2765 30152 43265 4132 3740 1210 1023 5935 1114
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
