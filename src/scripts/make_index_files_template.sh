#!/bin/bash

PIX_FILE=~mbusha/projects/development/final_healpix_cells.txt

cd 2MASS
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_2mass_mag."$1".fit\">DUMMY_NAME_2mass_mag."$1".fit</a><br>"}' > index.html
cd ..

cd BCS
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_bcs_mag."$1".fit\">DUMMY_NAME_bcs_mag."$1".fit</a><br>"}' > index.html
cd ..

cd CFHTLS
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_cfhtls_mag."$1".fit\">DUMMY_NAME_cfhtls_mag."$1".fit</a><br>"}' > index.html
cd ..

cd DEEP2
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_deep_mag."$1".fit\">DUMMY_NAME_deep_mag."$1".fit</a><br>"}' > index.html
cd ..

cd DR8
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_sdss_mag."$1".fit\">DUMMY_NAME_sdss_mag."$1".fit</a><br>"}' > index.html
cd ..

cd Euclid
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_euclid_mag."$1".fit\">DUMMY_NAME_euclid_mag."$1".fit</a><br>"}' > index.html
cd ..

cd FLAMEX
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_flamex_mag."$1".fit\">DUMMY_NAME_flamex_mag."$1".fit</a><br>"}' > index.html
cd ..

cd HSC
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_hsc_mag."$1".fit\">DUMMY_NAME_hsc_mag."$1".fit</a><br>"}' > index.html
cd ..

cd Johnson
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_johnson_mag."$1".fit\">DUMMY_NAME_johnson_mag."$1".fit</a><br>"}' > index.html
cd ..

cd LSST
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_lsst_mag."$1".fit\">DUMMY_NAME_lsst_mag."$1".fit</a><br>"}' > index.html
cd ..

cd NDWFS
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_ndwfs_mag."$1".fit\">DUMMY_NAME_ndwfs_mag."$1".fit</a><br>"}' > index.html
cd ..

cd RCS
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_rcs_mag."$1".fit\">DUMMY_NAME_rcs_mag."$1".fit</a><br>"}' > index.html
cd ..

cd Stripe82
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_stripe82_mag."$1".fit\">DUMMY_NAME_stripe82_mag."$1".fit</a><br>"}' > index.html
cd ..

cd VHS
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_vhs_mag."$1".fit\">DUMMY_NAME_vhs_mag."$1".fit</a><br>"}' > index.html
cd ..

cd VIKING
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_viking_mag."$1".fit\">DUMMY_NAME_viking_mag."$1".fit</a><br>"}' > index.html
cd ..

cd WFIRST
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_wfirst_mag."$1".fit\">DUMMY_NAME_wfirst_mag."$1".fit</a><br>"}' > index.html
cd ..

cd WISE
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_wise_mag."$1".fit\">DUMMY_NAME_wise_mag."$1".fit</a><br>"}' > index.html
cd ..

cd mask
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_mask."$1".fit\">DUMMY_NAME_mask."$1".fit</a><br>"}' > index.html
cd ..

cd photoz
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_zcarlos."$1".fit\">DUMMY_NAME_zCarlos."$1".fit</a><br>"}' > index.html
cd ..

cd photoz_DR8
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_DR8_zcarlos."$1".fit\">DUMMY_NAME_DR8_zCarlos."$1".fit</a><br>"}' > index.html
cd ..

cd ArborZ
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_arborz."$1".fit\">DUMMY_NAME_ArborZ."$1".fit</a><br>"}' > index.html
cd ..

cd obs_rotated
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME."$1".fit\">DUMMY_NAME."$1".fit</a><br>"}' > index.html
cd ..

cd truth_rotated
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_truth."$1".fit\">DUMMY_NAME_truth_des."$1".fit</a><br>"}' > index.html
cd ..

STARS_DIR=/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Stars/reformatted/
STARS_TBASE=${STARS_DIR}/Aardvark_0.5c_truth_stars
STARS_OBASE=${STARS_DIR}/Aardvark_0.5c_stars
mkdir -p Stars_truth
cd Stars_truth
cat $PIX_FILE | awk '{print "ln -s /nfs/slac/g/ki/ki19/des/mbusha/catalogs/Stars/reformatted//Aardvark_0.5c_truth_stars."$1".fit DUMMY_NAME_truth_stars."$1".fit"}' | sh
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_truth_stars."$1".fit\">DUMMY_NAME_truth_stars."$1".fit</a><br>"}' > index.html
cd ..

mkdir -p Stars_obs
cd Stars_obs
cat $PIX_FILE | awk '{print "ln -s /nfs/slac/g/ki/ki19/des/mbusha/catalogs/Stars/reformatted//Aardvark_0.5c_stars."$1".fit DUMMY_NAME_stars."$1".fit"}' | sh
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_stars."$1".fit\">DUMMY_NAME_stars."$1".fit</a><br>"}' > index.html
cd ..

mkdir -p ADDQSO_truth
cd ADDQSO_truth
cat $PIX_FILE | awk '{print "ln -s /nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/Brazil/v1.0/truth/ADDQSO_v1.0_truth."$1".fit DUMMY_NAME_addqso_truth."$1".fit"}' | sh
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_addqso_truth."$1".fit\">DUMMY_NAME_adqso_truth."$1".fit</a><br>"}' > index.html
cd ..

mkdir -p ADDQSO_obs
cd ADDQSO_obs
cat $PIX_FILE | awk '{print "ln -s /nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/Brazil/v1.0/obs/ADDQSO_v1.0."$1".fit DUMMY_NAME_addqso."$1".fit"}' | sh
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_addqso."$1".fit\">DUMMY_NAME_adqso."$1".fit</a><br>"}' > index.html
cd ..

mkdir -p ADDQSO_mask
cd ADDQSO_mask
cat $PIX_FILE | awk '{print "ln -s /nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/Brazil/v1.0/mask/ADDQSO_v1.0_mask."$1".fit DUMMY_NAME_addqso_mask."$1".fit"}' | sh
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_addqso_mask."$1".fit\">DUMMY_NAME_adqso_mask."$1".fit</a><br>"}' > index.html
cd ..

mkdir -p DESQSO_truth
cd DESQSO_truth
cat $PIX_FILE | awk '{print "ln -s /nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/DESQSO/truth/DESQSO_truth."$1".fit DUMMY_NAME_desqso_truth."$1".fit"}' | sh
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_desqso_truth."$1".fit\">DUMMY_NAME_desqso_truth."$1".fit</a><br>"}' > index.html
cd ..

mkdir -p DESQSO_obs
cd DESQSO_obs
cat $PIX_FILE | awk '{print "ln -s /nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/DESQSO/obs/DESQSO."$1".fit DUMMY_NAME_desqso."$1".fit"}' | sh
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_desqso."$1".fit\">DUMMY_NAME_desqso."$1".fit</a><br>"}' > index.html
cd ..

mkdir -p DESQSO_mask
cd DESQSO_mask
cat $PIX_FILE | awk '{print "ln -s /nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/DESQSO/mask/DESQSO_mask."$1".fit DUMMY_NAME_desqso_mask."$1".fit"}' | sh
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_desqso_mask."$1".fit\">DUMMY_NAME_desqso_mask."$1".fit</a><br>"}' > index.html
cd ..

mkdir -p DESQSO_VHS
cd DESQSO_VHS
cat $PIX_FILE | awk '{print "ln -s /nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/DESQSO/VHS/DESQSO."$1".fit DUMMY_NAME_desqso_vhs_mag."$1".fit"}' | sh
cat $PIX_FILE | awk '{print "<a href=\"DUMMY_NAME_desqso_vhs_mag."$1".fit\">DUMMY_NAME_desqso_vhs_mag."$1".fit</a><br>"}' > index.html
cd ..


cd halos
echo '<a href="DUMMY_NAME_halos_rotated.0.fit">DUMMY_NAME_halos.0.fit</a><br>' > index.html
echo '<a href="DUMMY_NAME_halos_rotated.1.fit">DUMMY_NAME_halos.1.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.4.fit">DUMMY_NAME_halos.4.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.5.fit">DUMMY_NAME_halos.5.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.6.fit">DUMMY_NAME_halos.6.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.7.fit">DUMMY_NAME_halos.7.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.12.fit">DUMMY_NAME_halos.12.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.13.fit">DUMMY_NAME_halos.13.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.14.fit">DUMMY_NAME_halos.14.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.15.fit">DUMMY_NAME_halos.15.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.16.fit">DUMMY_NAME_halos.16.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.20.fit">DUMMY_NAME_halos.20.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.21.fit">DUMMY_NAME_halos.21.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.22.fit">DUMMY_NAME_halos.22.fit</a><br>' >> index.html
echo '<a href="DUMMY_NAME_halos_rotated.23.fit">DUMMY_NAME_halos.23.fit</a><br>' >> index.html
cd ..

