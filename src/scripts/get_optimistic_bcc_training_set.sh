NAME=get_bcc_training_set
OUT=get_optimistic_bcc_training_set.out
ERR=get_optimistic_bcc_training_set.err
bsub -q kipac-ibq -oo $OUT -eo $ERR -J $NAME /afs/slac.stanford.edu/g/ki/software/idl/idl_6.3/bin/idl get_optimistic_bcc_training_set.idl

