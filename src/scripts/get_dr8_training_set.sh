NAME=get_dr8_training_set
OUT=$NAME.out
bsub -q kipac-ibq -oo $OUT -J $NAME /afs/slac.stanford.edu/g/ki/software/idl/idl_6.3/bin/idl get_dr8_training_set.idl
