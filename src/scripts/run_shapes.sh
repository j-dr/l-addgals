#!/bin/sh
BOX=$1
PIX=$2
NAME=lensing_${BOX}_${PIX}
OUT=shapes_logs/$NAME.out
ERR=shapes_logs/$NAME.err
bsub -q kipac-ibq -oo $OUT -eo $ERR -J $NAME ./make_shapes.sh $BOX $PIX
