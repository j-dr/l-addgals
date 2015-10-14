PIXNUM=$1
ZNUM=$2

cd $PIXNUM
cd ${ZNUM}
NAME=
OUT=${NAME}.out
ERR=${NAME}.err
rm -f $OUT
rm -f $ERR
echo $NUM
pwd
bsub -c 24:00 -a OpenMPI -q kipac-ibq -J $NAME -o $OUT -e $ERR -n 1 ../../run.sh ${NUM}
