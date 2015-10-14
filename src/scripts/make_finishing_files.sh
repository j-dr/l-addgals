#!/bin/bash

#variables that get passed
OUTPATH=$1
NAME=$2
SHEARBASE=$3
FINAL_NAME=$4

SCRIPTS_DIR=${OUTPATH}/scripts
mkdir -p $SCRIPTS_DIR

#the shapes file
SHAPES_TEMPLATE=make_shapes_template.sh
SHAPES_EXE=${SCRIPTS_DIR}/make_shapes.sh
RUN_SHAPES_TEMPLATE=run_shapes.sh
RUN_SHAPES_EXE=${SCRIPTS_DIR}/run_shapes.sh
RUN_ALL_SHAPES_TEMPLATE=run_all_shapes.sh
RUN_ALL_SHAPES_EXE=${SCRIPTS_DIR}/run_all_shapes.sh

sed -e 's:DUMMY_NAME:'PO_$NAME':'\
    < $SHAPES_TEMPLATE > $SHAPES_EXE
chmod 744 $SHAPES_EXE
cp $RUN_SHAPES_TEMPLATE $RUN_SHAPES_EXE
cp $RUN_ALL_SHAPES_TEMPLATE $RUN_ALL_SHAPES_EXE
mkdir -p ${SCRIPTS_DIR}/shapes_logs

#the finalize routines
FINALIZE_FILE=finalize_bcc_catalog.sav
RUN_FINALIZE_TEMPLATE=finalize_bcc_catalog_template.sh
RUN_FINALIZE=${SCRIPTS_DIR}/finalize_bcc_catalog.sh
RUN_ALL_FINALIZE=finalize_all_bcc_catalog.sh

mkdir -p ${SCRIPTS_DIR}/finalize_logs
cp $FINALIZE_FILE $SCRIPTS_DIR
cp $RUN_ALL_FINALIZE $SCRIPTS_DIR
sed -e 's:INPATH=:INPATH='$OUTPATH':'\
    -e 's:OUTPATH=:OUTPATH='$OUTPATH':'\
    -e 's:SHEARBASE=:SHEARBASE='$SHEARBASE':'\
    -e 's:ONAME=:ONAME='$FINAL_NAME':'\
    -e 's:FBASE=:FBASE=PO_'$NAME':'\
    < $RUN_FINALIZE_TEMPLATE > ${RUN_FINALIZE}
chmod 744 $RUN_FINALIZE

#the DR8 training set
DR8_TRAINING_SET=${OUTPATH}/photoz_DR8/${FINAL_NAME}_DR8_training_set.fit
DR8_TRAINING_SET2=${OUTPATH}/photoz_DR8/${FINAL_NAME}_DR8_training_set_sdss_mag.fit
DR8_TRAINING_TEMPLATE=get_dr8_training_set_template.idl
DR8_TRAINING=${SCRIPTS_DIR}/get_dr8_training_set.idl
RUN_DR8=get_dr8_training_set.sh

sed -e 's:OUTPATH=:outpath="'$OUTPATH'/photoz_DR8/":'\
    -e 's:PATH=:path="'$OUTPATH'":'\
    -e 's:TRUTH_BASE=:truth_base="'$FINAL_NAME'_truth":'\
    -e 's:SDSS_BASE=:sdss_base="'$FINAL_NAME'_sdss_mag":'\
    -e 's:OUTFILE1=:outfile1="'$DR8_TRAINING_SET':'\
    -e 's:OUTFILE2=:outfile2="'$DR8_TRAINING_SET2':'\
    < $DR8_TRAINING_TEMPLATE > $DR8_TRAINING
cp $RUN_DR8 ${SCRIPTS_DIR}

#the DES training set
DES_TRAINING_SET=${OUTPATH}/photoz/${FINAL_NAME}_Optimistic_training_set.fit
DES_TRAINING_TEMPLATE=get_optimistic_bcc_training_set_template.idl
DES_TRAINING=${SCRIPTS_DIR}/get_optimistic_bcc_training_set.idl
RUN_DES_TRAINING=get_optimistic_bcc_training_set.sh

sed -e 's:FBASE=:fbase = "'$FINAL_NAME'_truth":'\
    -e 's:OUTFILE=:outfile="'$DES_TRAINING_SET':'\
    < $DES_TRAINING_TEMPLATE > $DES_TRAINING
cp $RUN_DES_TRAINING ${SCRIPTS_DIR}

#the photo-z run/submission scripts
PHOTOZ_SAV_FILE=run_zcarlos_bcc.sav
PHOTOZ_RUN_TEMPLATE=run_zcarlos_bcc_template.sh
PHOTOZ_RUN_FILE=${SCRIPTS_DIR}/run_zcarlos_bcc.sh
PHOTOZ_RUN_ALL_FILE=run_all_zcarlos_bcc.sh
GFILE=${OUTPATH}/truth/${FINAL_NAME}_truth
DR8_PHOTOZDIR=${OUTPATH}/photoz_DR8
DES_PHOTOZDIR=${OUTPATH}/photoz
DR8_OBSFILE=${OUTPATH}/DR8/${FINAL_NAME}_sdss_mag
DR8_OFILE=${DR8_PHOTOZDIR}/${FINAL_NAME}_DR8_zcarlos
DES_OBSFILE=${OUTPATH}/truth/${FINAL_NAME}_truth
DES_OFILE=${DES_PHOTOZDIR}/${FINAL_NAME}_zcarlos
CHECK_DR8_PHOTOZ_FILE=check_dr8_photoz.sh
CHECK_PHOTOZ_FILE=check_photoz.sh

mkdir -p $DR8_PHOTOZDIR
mkdir -p $DES_PHOTOZDIR
mkdir -p ${SCRIPTS_DIR}/photoz_logs
cp $PHOTOZ_SAV_FILE $SCRIPTS_DIR
cp $PHOTOZ_RUN_ALL_FILE $SCRIPTS_DIR
cp $CHECK_DR8_PHOTOZ_FILE $SCRIPTS_DIR
cp $CHECK_PHOTOZ_FILE $SCRIPTS_DIR
sed -e 's:GFILE=:GFILE='$GFILE'.${PIXEL}.fit:'\
    -e 's:DR8_TRAINING_SET=:DR8_TRAINING_SET='$DR8_TRAINING_SET2':'\
    -e 's:DR8_OBSFILE=:DR8_OBSFILE='$DR8_OBSFILE'.${PIXEL}.fit:'\
    -e 's:DR8_OFILE=:DR8_OFILE='$DR8_OFILE'.${PIXEL}.fit:'\
    -e 's:DES_TRAINING_SET=:DES_TRAINING_SET='$DES_TRAINING_SET':'\
    -e 's:DES_OBSFILE=:DES_OBSFILE='$DES_OBSFILE'.${PIXEL}.fit:'\
    -e 's:DES_OFILE=:DES_OFILE='$DES_OFILE'.${PIXEL}.fit:'\
    < $PHOTOZ_RUN_TEMPLATE > $PHOTOZ_RUN_FILE
chmod 744 $PHOTOZ_RUN_FILE

#the rotation scripts
ROT_CAT=rotate_catalog.py
ROT_TOOLS=rot_mock_tools.py
ROT_SH=rotate_catalog.sh
ROT_ALL_TEMPLATE=rotate_all_catalogs_template.sh
ROT_ALL=${SCRIPTS_DIR}/rotate_all_catalogs.sh
ROT_OBSDIR=${OUTPATH}/obs_rotated
ROT_TRUTHDIR=${OUTPATH}/truth_rotated
RTIN=${OUTPATH}/truth/${FINAL_NAME}_truth
RTOUT=${ROT_TRUTHDIR}/${FINAL_NAME}_truth
ROIN=${OUTPATH}/obs/${FINAL_NAME}
ROOUT=${ROT_OBSDIR}/${FINAL_NAME}

mkdir -p $ROT_OBSDIR
mkdir -p $ROT_TRUTHDIR
mkdir -p ${SCRIPTS_DIR}/rotate_logs
cp $ROT_CAT $SCRIPTS_DIR
cp $ROT_TOOLS $SCRIPTS_DIR
cp $ROT_SH $SCRIPTS_DIR
sed -e 's:FIN_TRUTH=:FIN_TRUTH='$RTIN'.$i:'\
    -e 's:FOUT_TRUTH=:FOUT_TRUTH='$RTOUT'.$i:'\
    -e 's:FIN_OBS=:FIN_OBS='$ROIN'.$i:'\
    -e 's:FOUT_OBS=:FOUT_OBS='$ROOUT'.$i:'\
    < $ROT_ALL_TEMPLATE > $ROT_ALL
chmod 744 $ROT_ALL

#the mask scripts
MASK_SAV=mask_bcc_pixel.sav
RUN_MASK_TEMPLATE=mask_bcc_pixel_template.sh
RUN_MASK=${SCRIPTS_DIR}/mask_bcc_pixel.sh
MASK_ALL=mask_all_pixels.sh

mkdir -p ${SCRIPTS_DIR}/mask_logs
mkdir -p ${OUTPATH}/mask
cp $MASK_SAV ${SCRIPTS_DIR}
cp $MASK_ALL $SCRIPTS_DIR
sed -e 's:INBASE=:INBASE='$ROOUT':'\
    -e 's:OUTBASE=:OUBASE=../mask/'$FINAL_NAME':'\
    < $RUN_MASK_TEMPLATE > $RUN_MASK
chmod 744 $RUN_MASK

#single rotate and mask file
ROT_MASK_TEMPLATE=rotate_and_mask_template.sh
ROT_MASK=${SCRIPTS_DIR}/rotate_and_mask.sh
RUN_ROT_MASK=run_rotate_and_mask.sh
RUN_ALL_ROT_MASK=run_all_rotate_and_mask.sh
CHECK_MASK=check_mask.sh

mkdir -p ${SCRIPTS_DIR}/rotate_and_mask_logs
sed -e 's:FIN_TRUTH=:FIN_TRUTH='$RTIN':'\
    -e 's:FOUT_TRUTH=:FOUT_TRUTH='$RTOUT':'\
    -e 's:FIN_OBS=:FIN_OBS='$ROIN':'\
    -e 's:FOUT_OBS=:FOUT_OBS='$ROOUT':'\
    -e 's:MASK_OUT=:MASK_OUT=../mask/'${FINAL_NAME}_mask':'\
    < $ROT_MASK_TEMPLATE > $ROT_MASK
chmod 744 $ROT_MASK
cp $RUN_ROT_MASK $SCRIPTS_DIR
cp $RUN_ALL_ROT_MASK $SCRIPTS_DIR
cp $CHECK_MASK $SCRIPTS_DIR

#script for creating the index files
INDEX_TEMPLATE=make_index_files_template.sh
INDEX_SCRIPT=${OUTPATH}/make_index_files.sh
STARS_DIR=/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Stars/reformatted/
STARS_TRUTH=${STARS_DIR}/Aardvark_0.5c_truth_stars
STARS_OBS=${STARS_DIR}/Aardvark_0.5c_stars
ADDQSO_DIR=/nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/Brazil/v1.0/
ADDQSO_TRUTH=${ADDQSO_DIR}/truth/ADDQSO_v1.0_truth
ADDQSO_OBS=${ADDQSO_DIR}/obs/ADDQSO_v1.0
ADDQSO_MASK=${ADDQSO_DIR}/mask/ADDQSO_v1.0_mask
DESQSO_DIR=/nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/DESQSO/
DESQSO_TRUTH=${DESQSO_DIR}/truth/DESQSO_truth
DESQSO_OBS=${DESQSO_DIR}/obs/DESQSO
DESQSO_VHS=${DESQSO_DIR}/VHS/DESQSO
DESQSO_MASK=${DESQSO_DIR}/mask/DESQSO_mask

sed -e 's:DUMMY_NAME:'$FINAL_NAME':g'\
    -e 's:DUMMY_STARS_TRUTH:'$STARS_TRUTH':'\
    -e 's:DUMMY_STARS_OBS:'$STARS_OBS':'\
    -e 's:DUMMY_ADDQSO_TRUTH:'$ADDQSO_TRUTH':'\
    -e 's:DUMMY_ADDQSO_OBS:'$ADDQSO_OBS':'\
    -e 's:DUMMY_ADDQSO_MASK:'$ADDQSO_MASK':'\
    -e 's:DUMMY_DESQSO_TRUTH:'$DESQSO_TRUTH':'\
    -e 's:DUMMY_DESQSO_OBS:'$DESQSO_OBS':'\
    -e 's:DUMMY_DESQSO_VHS:'$DESQSO_VHS':'\
    -e 's:DUMMY_DESQSO_MASK:'$DESQSO_MASK':'\
    < $INDEX_TEMPLATE > $INDEX_SCRIPT
chmod 744 $INDEX_SCRIPT

#script for creating the halo catalog
HALO_SCRIPT=make_halo_catalog.sh
HALO_IDL_TEMPLATE=make_halo_catalog_template.idl
HALO_IDL=${SCRIPTS_DIR}/make_halo_catalog.idl

mkdir -p ${OUTPATH}/halos
cp $HALO_SCRIPT $SCRIPTS_DIR
sed -e 's:PATH=:path = "'$OUTPATH'/individual_box_files/":'\
    -e 's:OUTNAME1=:outname1 = "PO_'$NAME'_1050_halos":'\
    -e 's:OUTNAME2=:outname2 = "PO_'$NAME'_2600_halos":'\
    -e 's:OUTNAME3=:outname3 = "PO_'$NAME'_4000_halos":'\
    -e 's:OUTBASE=:outbase = "../halos/'$FINAL_NAME'_halos":'\
    < $HALO_IDL_TEMPLATE > $HALO_IDL
