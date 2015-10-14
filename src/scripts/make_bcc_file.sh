PIX=$1
FPATH=
FBASE=
OUTPATH=
HFILE=
LABEL=
NFILES=

NAME=make_bcc_file_${LABEL}.$PIX
OUT=make_bcc_file_logs/$NAME.out
ERR=make_bcc_file_logs/$NAME.err
bsub -n 4 -q kipac-ibq -c 24:00 -oo $OUT -eo $ERR -J $NAME idl_vm_run.py -s make_bcc_file.sav $PIX $FPATH $FBASE $OUTPATH $HFILE $NFILES
