#!/bin/bash

PIX=$1
NAME=finalize.$PIX
OUT=finalize_logs/$NAME.out
ERR=finalize_logs/$NAME.err

INPATH=
OUTPATH=
SHEARBASE=
DO_LENSING=1
ONAME=
FBASE=

bsub -x -q kipac-ibq -oo $OUT -eo $ERR -J $NAME idl_vm_run.py finalize_bcc_catalog.sav $INPATH $OUTPATH $ONAME $SHEARBASE $PIX $DO_LENSING $FBASE
