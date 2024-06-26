#! /bin/bash

SOURCE=$1
MATING_FUNCTION=$2
OUTDIR=$3

LEAD='^## BEGIN GENERATED CONTENT$'
TAIL='^## END GENERATED CONTENT$'

SLIM_MODEL=$(sed -e "/$LEAD/,/$TAIL/{ /$LEAD/{r $MATING_FUNCTION
  d;}; /$TAIL/d; }" $SOURCE)

echo "$SLIM_MODEL" > $OUTDIR/model.slim 
