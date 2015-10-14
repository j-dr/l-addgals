nboxes = 4
sim_zmin = [0.0, 0.32, 0.87, 1.97]
sim_zmax = [0.36, 0.93, 2.04, 5.0]
bcg_mass_lim = ['1.35e12', '6.5e12', '1e13', '4.0e13']
nproc = [8,20,18, 20]
omegam = 0.23
omegal = 0.77
boxsize = ['1050', '2600', '4000', '6000']
simname = 'Aardvark'
;simfile = '/lustre/ki/pfs/mbusha/projects/Aardvark/Lb'+boxsize+$
;	   '/output/snapshot_Lightcone_000_healpix'
simfile = '/nfs/slac/g/ki/ki21/cosmo/jderose/simulations/lightcone/Aardvark/Lb'+boxsize+$
	   '/snapshot_Lightcone_000_healpix'
rnnfile = '/nfs/slac/g/ki/ki21/cosmo/jderose/simulations/rnn/Aardvark/Lb'+boxsize+$
	   '/rnn_snapshot_Lightcone_000_healpix'
;rnnfile = '/lustre/ki/pfs/mbusha/projects/Aardvark/Lb'+boxsize+$
;	  '/rnn/rnn_snapshot_Lightcone_000_healpix'
halofile = '/nfs/slac/g/ki/ki21/cosmo/jderose/simulations/halos/Aardvark/Lb'+boxsize+$
	   '/out_0.parents'
;halofile = '/lustre/ki/pfs/mbusha/projects/Aardvark/Lb'+boxsize+$
;	   '/Quartant1/rockstar/out_0.parents'
hfile = '/nfs/slac/g/ki/ki21/cosmo/jderose/simulations/halos/Aardvark/Lb'+boxsize+$
	   '/halos'+boxsize+'.fit'
;hfile = '/lustre/ki/pfs/mbusha/projects/Aardvark/Lb'+boxsize+$
;	   '/Quartant1/rockstar/halos_'+boxsize+'.fit'
;rnn_halofile = '/nfs/slac/g/ki/ki19/des/mbusha/simulations/Aardvark-2PS/Lb'+boxsize+'/Quartant1/calcrnn_halos/outputs/rnn_out_0.parents'
;rnn_halofile =
;'/lustre/ki/pfs/mbusha/projects/Aardvark/Lb'+boxsize+'/Quartant1/calcrnn_halos/outputs/rnn_out_0.parents'
rnn_halofile = '/nfs/slac/g/ki/ki21/cosmo/jderose/simulations/halos/Aardvark/Lb'+boxsize+$
           '/rnn_out_0.parents'
dir = '/nfs/slac/g/ki/ki21/cosmo/jderose/addgals/catalogs/Aardvark/Lb'+boxsize+'_v1.0'
catdir = '/nfs/slac/g/ki/ki21/cosmo/jderose/addgals/catalogs/Aardvark/Catalog_v1.0/'
;;;I don't think I'm still using denspdfstr
;denspdfstr = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchilla-tuning/analysis/rdel/denspdf_Chinchilla.dat'
;denspdfstr=['/u/ki/mbusha/projects/modules/idl/addgals/rdel/denspdf_Consuelo02_AGES_rescale_099.dat', 

;for i = 0, nboxes - 1 do begin
;for i = 0, nboxes -  2  do begin
i=2
make_bcc_box_aardvark, dir=dir[i], $
                       sim_zmin=sim_zmin[i], sim_zmax=sim_zmax[i], $
                       nproc=nproc[i], $
                       omegam=omegam, omegal=omegal, $
                       simfile=simfile[i], rnnfile=rnnfile[i], $
                       halofile=halofile[i], rnn_halofile=rnn_halofile[i], $
                       simname=simname, boxsize=boxsize[i], $
                       hfile=hfile[i], bcg_mass_lim=bcg_mass_lim[i],catdir=catdir
;endfor

end
