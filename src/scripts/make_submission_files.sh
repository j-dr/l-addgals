#!/bin/bash

SIMNAME=$1
BOXSIZE=$2
DIR=$3

#our template file names
RUN=run.sh
RUNSCRIPT=runscript.sh
SUB=submit_jobs.sh
CHECK_GINFO=check_ginfo1_files.sh
RESUB_GINFO=check_and_finish_ginfo1_files.sh
CHECK_IDL=check_idl_files.sh
RESUB_IDL=check_and_finish_idl.sh

#the submission script
case $BOXSIZE in
	1050)
		ZSTR='000 001 002 003 004 005 006 007'
		NPROC=1
		;;
	2600)
		ZSTR='000 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019'
		NPROC=2
		;;
	4000)
		ZSTR='000 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017'
		NPROC=2
		;;
	6000)
		ZSTR='000 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019'
                NPROC=2
                ;;
	*)
esac

cp $RUN $DIR
sed -e 's:NAME=:NAME='$SIMNAME'_'$BOXSIZE'.${PIXNUM}.${ZNUM}:'\
    -e 's:-n 1:-n '$NPROC':'\
	<$RUNSCRIPT > $DIR/$RUNSCRIPT
chmod 744 $DIR/$RUNSCRIPT
sed -e 's:for Z in:for Z in '"$ZSTR"':'\
        <$SUB > $DIR/$SUB
chmod 744 $DIR/$SUB

#the files for checking completion
sed -e 's:for Z in:for Z in '"$ZSTR"':'\
        <$CHECK_GINFO > $DIR/$CHECK_GINFO
chmod 744 $DIR/$CHECK_GINFO
sed -e 's:for Z in:for Z in '"$ZSTR"':'\
        <$RESUB_GINFO > $DIR/$RESUB_GINFO
chmod 744 $DIR/$RESUB_GINFO
sed -e 's:for Z in:for Z in '"$ZSTR"':'\
    -e 's:PO_:PO_'$SIMNAME'_'$BOXSIZE':'\
	<$CHECK_IDL > $DIR/$CHECK_IDL	
chmod 744 $DIR/$CHECK_IDL
sed -e 's:for Z in:for Z in '"$ZSTR"':'\
    -e 's:PO_:PO_'$SIMNAME'_'$BOXSIZE':'\
        <$RESUB_IDL > $DIR/$RESUB_IDL
chmod 744 $DIR/$RESUB_IDL
