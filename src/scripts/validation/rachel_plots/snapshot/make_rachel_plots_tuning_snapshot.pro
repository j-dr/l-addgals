labelarr = ['bcgs_1e12', 'bcgs_1e13', 'bcgs_1e14', 'no_bcgs', 'bcgs_3e12', 'bcgs_3e13', 'bcgs_3e14']
;labelarr = ['bcgs_1e14', 'no_bcgs', 'bcgs_3e12', 'bcgs_3e13', 'bcgs_3e14']

for i = 1L, N_ELEMENTS(labelarr) - 1 do begin
  label = labelarr[i]

  dir = '/nfs/slac/g/ki/ki23/des/mbusha/Simulations/Chinchilla-tuning/downsample/addgals/single_epoch_fit/'+label+'/idl/'
  outdir = dir+'plot_info/'
  gfname = 'PO_Chinchilla-1_1050_0.000_galaxies.fit'
  hfname = 'PO_Chinchilla-1_1050_0.000_halos.fit'
  masscut = 1e12

  make_rachel_plotdata_snapshot, outdir, dir, gfname=gfname, hfname=hfname, $
        masscut=masscut
endfor

end
