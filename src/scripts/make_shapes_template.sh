#!/bin/sh
LB=$1
PIXEL=$2

TRUTH=DUMMY_NAME_${LB}_truth_no_photoz.${PIXEL}.fit

echo $TRUTH

SHAPES_EXE=~/projects/mockshapes/mockshapes-0.5.0/src/mockshapes
LENSING_EXE=~/projects/mockshapes/mockshapes-0.5.0/src/apply_lensing
IN_PATH=../shapes/
TPATH=../individual_box_files/
RANDOM=`date '+%s'`

###make copies of all the files
echo "Making copies..."
mkdir $IN_PATH
cp ${TPATH}/${TRUTH} $IN_PATH

###now add the shapes 
echo "Adding shapes..."
$SHAPES_EXE -SEED $RANDOM $IN_PATH/$TRUTH
echo "Succesfully finished mockshapes $LB $PIXEL!"

