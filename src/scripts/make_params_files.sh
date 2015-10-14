#!/bin/sh

ZMIN=0.1
ZMAX=0.2
NUM=10
MMAG=-18.0

ZMIN=$1
ZMAX=$2
OPATH=$3
MMAG=$4
DENSFILE=$5
LBCGFILE=$6
PHISTAR=$7

STRINGFILE=StringParameters
NUMFILE=NumericalParameters

cp StringParameters_maxbcg_template $STRINGFILE
cp NumericalParameters_maxbcg_template $NUMFILE

echo "ZREDMIN $ZMIN" >> $NUMFILE
echo "ZREDMAX $ZMAX" >> $NUMFILE
echo "Magmin $MMAG" >> $NUMFILE
echo "out_path $OPATH" >> $STRINGFILE 
echo "path $OPATH" >> $STRINGFILE 
echo "denspdffile $DENSFILE" >> $STRINGFILE
echo "lbcgfile $LBCGFILE" >> $STRINGFILE
echo "phistar $PHISTAR" >> $NUMFILE
