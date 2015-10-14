#!/bin/sh

FPATH=$1
FBASE=$2
OUTPATH=$3
HFILE=$4
NFILE=$5

# now we put together the file to combine the redshift bins
MAKE_BCC_SAV=make_bcc_file.sav
MAKE_BCC_SH=make_bcc_file.sh

echo " "

# copy the .sav file over
cp $MAKE_BCC_SAV ${FPATH}

# make the actual submission script
sed -e 's:FPATH=:FPATH='$FPATH':'\
    -e 's:FBASE=:FBASE='$FBASE':'\
    -e 's:OUTPATH=:OUTPATH='$OUTPATH':'\
    -e 's:HFILE=:HFILE='$HFILE':'\
    -e 's:LABEL=:LABEL='$FBASE':'\
    -e 's:NFILES=:NFILES='$NFILE':'\
  < $MAKE_BCC_SH > ${FPATH}/${MAKE_BCC_SH}
chmod 744 ${FPATH}/${MAKE_BCC_SH}
