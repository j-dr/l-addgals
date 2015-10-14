; the number of boxes we want to make
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
omegam = 0.286
omegal = 0.714

; the file containing the denspdf parameterization
paramfile = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla_v2/rdel/parameters'+['_Lb1050', '', '', '']+'.txt'

;;; files where all the information is sored
 ; simulation name
simname = 'Chinchilla-1'

 ; files for the particles and rnn information
simfile = '/lustre/ki/pfs/mbusha/projects/Chinchilla-1/Lb'+boxsize+$
	   '/output/snapshot_Lightcone_000_healpix'
rnnfile = '/lustre/ki/pfs/mbusha/projects/Chinchilla-1/Lb'+boxsize+$
	  '/rnn/rnn_snapshot_Lightcone_000_healpix'

 ; halo file -- one is a rockstar file, one is in .fit format
halofile = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb'+boxsize+'/Octants_01/rockstar/out_0.parents'
hfile = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb'+boxsize+'/Octants_01/rockstar/halos_'+boxsize+'.fit'
rnn_halofile = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb'+boxsize+'/Octants_01/calcrnn_halos/outputs/rnn_out_0.parents'

 ; the directory where the cosmological box is stored
dir = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla-1/Lb'+boxsize+'_v2.0'

 ; the directory where we store the final catalog
catdir = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla-1/Catalog_v2.0/'

;denspdfstr = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchilla-tuning/analysis/rdel/denspdf_Chinchilla.dat'
denspdfstr = '/kaka'

for i = 0, nboxes - 1 do begin
  make_bcc_box, dir=dir[i], $
        sim_zmin=sim_zmin[i], sim_zmax=sim_zmax[i], $
        nproc=nproc[i], $
        omegam=omegam, omegal=omegal, $
        simfile=simfile[i], rnnfile=rnnfile[i], $
        halofile=halofile[i], rnn_halofile=rnn_halofile[i], $
        simname=simname, boxsize=boxsize[i], denspdfstr=denspdfstr, $
	bcg_mass_lim = bcg_mass_lim[i], paramfile=paramfile[i], $
	catdir=catdir, hfile=hfile[i];, npix=1
endfor

end
