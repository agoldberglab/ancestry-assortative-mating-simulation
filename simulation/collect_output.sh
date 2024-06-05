#! /bin/bash

OUTPUT=/path/to/simulation/output
SUMMARY=$OUTPUT/summary

TABLES=(binned_ancestry0.csv
	shuffled_rho.csv
	shuffled_delta.csv
	summary.csv)

mkdir -p $SUMMARY

for TABLE in ${TABLES[@]}; do
    > $SUMMARY/$TABLE
    
    COUNTER=0
    
    while IFS=, read -r SEED MATE_CHOICE THETA INITIAL \
	MIG_RATE POPULATION_SIZE; do
    
	THETA="likelihood_"$THETA"x"
	OUTDIR=$THETA"/"$MATE_CHOICE"/migration_rate_"$MIG_RATE"/initial_"$INITIAL"/"$SEED

        FILENAME=$OUTPUT/$OUTDIR/$TABLE
        COUNTER=$((COUNTER+1))

        if [ $COUNTER == 1 ]; then
	    cat $FILENAME > tmp.csv
	    awk -v s=$COUNTER 'BEGIN { FS=OFS="," } \
            { if (NR==1) { print "simulation" OFS $0 } \
              else { print s OFS $0} } ' tmp.csv >> $SUMMARY/$TABLE
        else
	    tail -n +2 $FILENAME > tmp.csv
            awk -v s=$COUNTER 'BEGIN { FS=OFS="," } \
              { print s OFS $0} ' tmp.csv >> $SUMMARY/$TABLE
	fi
    
    done < parameters.csv #SIMULATIONS

done #TABLES
rm tmp.csv
