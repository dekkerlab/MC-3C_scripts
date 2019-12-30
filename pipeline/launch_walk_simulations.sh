

for i in $(seq 17122016 1 17122116)
do
	bsub -q long -W 480 -n 24 -R "rusage[mem=64000]" "Rscript mitotic_simulation_test.R $i; Rscript ns_simulation_test.R $i; Rscript interactions_to_usable_frame_for_sims.R $i"
done

