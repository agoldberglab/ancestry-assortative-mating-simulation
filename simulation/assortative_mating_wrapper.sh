#! /bin/bash

SCRIPTS=/path/to/simulation/scripts
OUTPUT=/path/to/simulation/output

> commands.txt

while IFS=, read -r SEED MATE_CHOICE THETA INITIAL \
	MIGRATION_RATE POPULATION_SIZE N_GENERATIONS; do

    LIKELIHOOD="likelihood_"$THETA"x"
    MIG_RATE="migration_rate_"$MIGRATION_RATE
    ADMIXTURE="initial_"$INITIAL
    
    OUTDIR=$LIKELIHOOD/$MATE_CHOICE/$MIG_RATE/$ADMIXTURE/$SEED

    CMD=$(echo python3 $SCRIPTS/assortative_mating_driver.py \
      --source $SCRIPTS/global_ancestry_assortative_mating.slim \
      --seed $SEED \
      --theta $THETA \
      --mate_choice $MATE_CHOICE \
      --recombination_rate_file \
          $SCRIPTS/twenty_two_chromosome_uniform_recombination_map.txt \
      --initial_proportion $INITIAL \
      --migration_rate $MIGRATION_RATE \
      --population_size $POPULATION_SIZE \
      --n_generations $N_GENERATIONS \
      --outdir $OUTPUT"/"$OUTDIR)
    
    echo $CMD >> commands.txt

done < parameters.csv
