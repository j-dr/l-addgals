nboxes = 3
sim_zmin = [0.0, 0.32, 0.87]
sim_zmax = [0.36, 0.93, 2.04]
nproc = [8,20,18]
omegam = 0.26
omegal = 0.74
boxsize = ['1050', '2600', '4000']
simname = 'Dingo'
simfile = '/lustre/ki/pfs/mbusha/projects/Dingo/Lb'+boxsize+$
	   '/output/snapshot_Lightcone_000_healpix'
rnnfile = '/lustre/ki/pfs/mbusha/projects/Dingo/Lb'+boxsize+$
	  '/rnn/rnn_snapshot_Lightcone_000_healpix'
halofile = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Dingo/Lb'+boxsize+'/Octants_01/rockstar/out_0.parents'
rnn_halofile = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Dingo/Lb'+boxsize+'/Octants_01/calcrnn_halos/outputs/rnn_out_0.parents'
dir = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Dingo/Lb'+boxsize+'_v1.0'

for i = 0, nboxes - 1 do begin
  make_bcc_box, dir=dir[i], $
        sim_zmin=sim_zmin[i], sim_zmax=sim_zmax[i], $
        nproc=nproc[i], $
        omegam=omegam, omegal=omegal, $
        simfile=simfile[i], rnnfile=rnnfile[i], $
        halofile=halofile[i], rnn_halofile=rnn_halofile[i], $
        simname=simname, boxsize=boxsize[i]
endfor

end
