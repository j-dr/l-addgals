; the number of boxes we want to make
pro make_fox
nboxes = 3

; the redshfit ranges for each of the boxes
sim_zmin = [0.0, 0.32, 0.87, 1.97]
sim_zmax = [0.36, 0.93, 2.04, 5.0]

; basic box properties -- these tend not to change
boxsize = ['1050', '2600', '4000', '6000']
bcg_mass_lim = ['1.35e12', '6.5e12', '1e13', '4.0e13']

; the number of redshift bins in each box
nproc = [8,20,18,20]

; cosmology
omegam = 0.31834
omegal = 0.68166

; the file containing the denspdf parameterization
paramfile = '/nfs/slac/g/ki/ki21/cosmo/jderose/halos/sham_and_histograms/catalogs/Fox_v2/rdel/parameters'+['_Lb1050', '_Lb2600', '_Lb4000', '_Lb6000']+'.txt'

;;; files where all the information is sored
 ; simulation name
simname = 'Fox-1'

 ; files for the particles and rnn information
simfile = '/nfs/slac/g/ki/ki22/cosmo/jderose/halos/simulation_preprocessing/FLb'+boxsize+$
	   '_output/snapshot_Lightcone_000_healpix'
rnnfile = '/nfs/slac/g/ki/ki22/cosmo/jderose/halos/simulation_preprocessing/FLb'+boxsize+$
	   '_rnn/snapshot_Lightcone_000_healpix'

 ; halo file -- one is a rockstar file, one is in .fit format
halofile = '/nfs/slac/g/ki/ki22/cosmo/jderose/skycats/scripts/rockstar/output/Lb'+boxsize+'/out_0.parents'
hfile = '/nfs/slac/g/ki/ki21/cosmo/jderose/halos/simulation_preprocessing/Fox/Lb'+boxsize+'/out_0.fits'
rnn_halofile = '/nfs/slac/g/ki/ki21/cosmo/jderose/halos/calcrnn/output/FLb'+boxsize+'/rnn_out_0.parents'

 ; the directory where the cosmological box is stored
dir = '/nfs/slac/g/ki/ki22/cosmo/jderose/addgals/catalogs/Fox/Lb'+boxsize

 ; the directory where we store the final catalog
catdir = '/nfs/slac/g/ki/ki22/cosmo/jderose/addgals/catalogs/Fox/Catalog_v1.0/'

denspdfstr = '/nfs/slac/g/ki/ki21/cosmo/jderose/halos/sham_and_histograms/catalogs/Fox_v2/rdel_for_addgals/parameters_Lb'+boxsize+'.txt'
;denspdfstr = '/kaka'

i=0
make_bcc_box, dir=dir[i], sim_zmin=sim_zmin[i], sim_zmax=sim_zmax[i], nproc=nproc[i], omegam=omegam, omegal=omegal, simfile=simfile[i], rnnfile=rnnfile[i], $
halofile=halofile[i], rnn_halofile=rnn_halofile[i], simname=simname, boxsize=boxsize[i], denspdfstr=denspdfstr, bcg_mass_lim = bcg_mass_lim[i], $
paramfile=paramfile[i], catdir=catdir, hfile=hfile[i]


end
