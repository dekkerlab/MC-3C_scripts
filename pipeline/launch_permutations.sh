
for i in $(seq 20161217 20161316); do bsub -R rusage[mem=64000] -n 32 -W 240 -q short "Rscript permutation_proof_of_principle.R $i; Rscript interactions_to_usable_frame_for_permutations.R $i"; done

