#!/bin/sh
PIX=$1
NAME=rotate_mask.${PIX}
OUT=rotate_and_mask_logs/$NAME.out
ERR=rotate_and_mask_logs/$NAME.err
bsub -q kipac-ibq -oo $OUT -eo $ERR -J $NAME rotate_and_mask.sh $PIX
