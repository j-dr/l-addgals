DUMMY_INBASE=/nfs/slac/g/ki/ki18/des/beckermr/aardvark_rotcat_v1.0_lowzfix/rotcats/obs/Aardvark_v0.5_des_masked
DUMMY_OUTBASE=../mask/Aardvark_v0.5_mask
ISTR=$1

echo "bsub -q kipac-ibq -oo mask_logs/mask_pixel_${ISTR}.out -J mask_pixel_$ISTR idl_vm_run.py mask_bcc_pixel.sav $INBASE $OUTBASE $ISTR"
