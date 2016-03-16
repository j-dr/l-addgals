#!/bin/sh

ZMIN=$1
ZMAX=$2
OPATH=$3
MMAG=$4
DENSFILE=$5
LBCGFILE=$6
PHISTAR=$7
PIXELNUM=$8
ZNUM=$9
OUTDIR=${10}
HALOFILE=${11}
RNN_HALOFILE=${12}
SIMNAME=${13}
BOXSIZE=${14}
PARAMFILE=${15}
BCG_MASS_LIM=${16}
DDIR=${17}
SRCDIR=${18}
FLABEL=${SIMNAME}_${BOXSIZE}

echo ${OUTDIR}
echo ${DDIR}
#echo $ZMIN, $ZMAX
#echo $OPATH
#echo $MMAG
#echo $DENSFILE
#echo $LBCGFILE
#echo $PHISTAR, $PIXELNUM, $ZNUM
#echo $SIMNAME
#echo $BOXSIZE

STRINGFILE=StringParameters
NUMFILE=NumericalParameters


cp $SRCDIR/scripts/StringParameters_BCC_template $STRINGFILE
cp $SRCDIR/scripts/NumericalParameters_BCC_template $NUMFILE
#cp StringParameters_BCC_template $STRINGFILE
#cp NumericalParameters_faber_4000_template.Buzzard $NUMFILE

# set the numberical parameters
echo "ZREDMIN $ZMIN" >> $NUMFILE
echo "ZREDMAX $ZMAX" >> $NUMFILE
echo "Magmin $MMAG" >> $NUMFILE
echo $PARAMFILE >> $NUMFILE
echo "BCG_Mass_lim $BCG_MASS_LIM" >> $NUMFILE

# set the string parameters
echo "datadir $DDIR" >> $STRINGFILE
echo "flabel $FLABEL" >> $STRINGFILE
echo "out_path $OPATH" >> $STRINGFILE 
echo "path $OPATH" >> $STRINGFILE 
echo "denspdffile $DENSFILE" >> $STRINGFILE
echo "lbcgfile $LBCGFILE" >> $STRINGFILE

echo "halofile  ${HALOFILE}" >> $STRINGFILE
echo "rnn_halofile         ${RNN_HALOFILE}" >> $STRINGFILE
echo "pstr ${PIXELNUM}" >> $STRINGFILE
echo "zstr ${ZNUM}" >> $STRINGFILE

echo "phistar $PHISTAR" >> $NUMFILE
echo "PixelNum $PIXELNUM" >> $NUMFILE

mkdir -p ${OUTDIR}
cp ${SRCDIR}/hv $OUTDIR
cp StringParameters $OUTDIR
cp NumericalParameters $OUTDIR

