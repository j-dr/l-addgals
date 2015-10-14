#!/bin/bash

###out input file names
PIX=$1
FIN_OBS=
FOUT_OBS=
FIN_TRUTH=
FOUT_TRUTH=
MASK_OUT=

#do the rotation
echo "Rotating the obs catalog..."
python rotate_catalog.py ${FIN_OBS}.${PIX}.fit ${FOUT_OBS}.${PIX}.fit
echo "Rotating the truth catalog..."
python rotate_catalog.py ${FIN_TRUTH}.${PIX}.fit ${FOUT_TRUTH}.${PIX}.fit

#do the masking
echo "Masking the obs catalog..."
idl_vm_run.py mask_bcc_pixel.sav $FOUT_OBS $MASK_OUT $PIX

echo "Completed rotate_and_mask for pixel " $PIX
