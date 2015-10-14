NAME=get_dr8_full_sky_training_set
OUT=$NAME.out
ERR=$NAME.err
bsub -q kipac-ibq -oo $OUT -eo $ERR -J $NAME /afs/slac.stanford.edu/g/ki/software/idl/idl_6.3/bin/idl get_dr8_full_sky_training_set.idl
