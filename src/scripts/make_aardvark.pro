nboxes = 4
sim_zmin = [0.0, 0.32, 0.87, 1.97]
sim_zmax = [0.36, 0.93, 2.04, 5.0]
nproc = [8,20,18, 20]
omegam = 0.23
omegal = 0.77
boxsize = ['1050', '2600', '4000', '6000']
simname = 'Aardvark'
simfile = '/lustre/ki/pfs/mbusha/projects/Aardvark/Lb'+boxsize+$
	   '/output/snapshot_Lightcone_000_healpix'
rnnfile = '/lustre/ki/pfs/mbusha/projects/Aardvark/Lb'+boxsize+$
	  '/rnn/rnn_snapshot_Lightcone_000_healpix'
halofile = '/lustre/ki/pfs/mbusha/projects/Aardvark/Lb'+boxsize+$
	   '/Quartant1/rockstar/out_0.parents'
;rnn_halofile = '/nfs/slac/g/ki/ki19/des/mbusha/simulations/Aardvark-2PS/Lb'+boxsize+'/Quartant1/calcrnn_halos/outputs/rnn_out_0.parents'
rnn_halofile = '/lustre/ki/pfs/mbusha/projects/Aardvark/Lb'+boxsize+'/Quartant1/calcrnn_halos/outputs/rnn_out_0.parents'
dir = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Aardvark/Lb'+boxsize+'_v1.01'
;;;I don't think I'm still using denspdfstr
denspdfstr = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchilla-tuning/analysis/rdel/denspdf_Chinchilla.dat'

;for i = 0, nboxes - 1 do begin
for i = 2, 2 do begin
  make_bcc_box, dir=dir[i], $
        sim_zmin=sim_zmin[i], sim_zmax=sim_zmax[i], $
        nproc=nproc[i], $
        omegam=omegam, omegal=omegal, $
        simfile=simfile[i], rnnfile=rnnfile[i], $
        halofile=halofile[i], rnn_halofile=rnn_halofile[i], $
        simname=simname, boxsize=boxsize[i],denspdfstr=denspdfstr;,npix=1
endfor

end
