PIXEL=$1
NAME=zCarlos_BCC.$PIXEL
OUT=../photoz_logs/${NAME}.out
GFILE=
DR8_TRAINING_SET=
DR8_OBSFILE=
DR8_OFILE=
DES_TRAINING_SET=
DES_OBSFILE=
DES_OFILE=

mkdir $PIXEL
cd $PIXEL
cp ../run_zcarlos_bcc.sav .

bsub -q kipac-ibq -oo $OUT -eo $OUT -J $NAME idl_vm_run.py ./run_zcarlos_bcc.sav $PIXEL $GFILE $DR8_TRAINING_SET $DR8_OBSFILE $DR8_OFILE $DES_TRAINING_SET $DES_OBSFILE $DES_OFILE

cd ..
