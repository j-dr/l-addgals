IN_NAME=$1
OUT_NAME=$2
NAME=$3
bsub -q kipac-ibq -oo rotate_logs/${NAME}.out -eo rotate_logs/${NAME}.err -J $NAME python rotate_catalog.py ${IN_NAME}.fit ${OUT_NAME}.fit
