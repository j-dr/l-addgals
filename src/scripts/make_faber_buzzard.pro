nboxes = 4
sim_zmin = [0.0, 0.32, 0.87, 1.97]
sim_zmax = [0.36, 0.93, 2.04, 6.0]
nproc = [8,20,18, 20]
omegam = 0.286
omegal = 0.714
boxsize = ['1050', '2600', '4000', '6000']
bcg_mass_lim = ['1.35e12', '6.5e12', '1e13', '4.0e13']
simname = 'Buzzard'
;simfile = '/lustre/ki/pfs/mbusha/projects/Chinchilla-1/Lb'+boxsize+$
;	   '/output/snapshot_Lightcone_000_healpix'
simfile = '/nfs/slac/g/ki/ki21/cosmo/jderose/simulations/lightcone/Buzzard/Lb'+boxsize+$
	   '/snapshot_Lightcone_000_healpix'
;rnnfile = '/lustre/ki/pfs/mbusha/projects/Chinchilla-1/Lb'+boxsize+$
;	  '/rnn/rnn_snapshot_Lightcone_000_healpix'
rnnfile = '/nfs/slac/g/ki/ki21/cosmo/jderose/simulations/rnn/Buzzard/Lb'+boxsize+$
	   '/rnn_snapshot_Lightcone_000_healpix'
halofile = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb'+boxsize+'/Octants_01/rockstar/out_0.parents'
hfile = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb'+boxsize+'/Octants_01/rockstar/halos_'+boxsize+'.fit'
rnn_halofile = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb'+boxsize+'/Octants_01/calcrnn_halos/outputs/rnn_out_0.parents'
dir = '/nfs/slac/g/ki/ki21/cosmo/jderose/addgals/catalogs/Buzzard/Lb'+boxsize+'_v1.1'
catdir = '/nfs/slac/g/ki/ki21/cosmo/jderose/addgals/catalogs/Buzzard/Catalog_v1.1/'
denspdfstr = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchilla-tuning/analysis/rdel/denspdf_Chinchilla.dat'

i=2
make_bcc_box_buzzard, dir=dir[i], $
              sim_zmin=sim_zmin[i], sim_zmax=sim_zmax[i], $
              nproc=nproc[i], $
              omegam=omegam, omegal=omegal, $
              simfile=simfile[i], rnnfile=rnnfile[i], $
              halofile=halofile[i], rnn_halofile=rnn_halofile[i], $
              simname=simname, boxsize=boxsize[i], $
              hfile=hfile[i], bcg_mass_lim=bcg_mass_lim[i],catdir=catdir

end
