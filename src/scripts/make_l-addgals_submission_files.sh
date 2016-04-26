#!/bin/bash

SIMNAME=$1
BOXSIZE=$2
DIR=$3


#our template file names
RUNSCRIPT=run_cell.sh
SUB=submit_jobs.sh
CHECK_GINFO=check_ginfo1_files.sh
RESUB_GINFO=check_and_finish_ginfo1_files.sh

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

sed -e 's:CATPATH=:CATPATH='$DIR':'\
    -e 's:>&:>&'$SIMNAME'_'$BOXSIZE'.${1}.${2}:'\
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
