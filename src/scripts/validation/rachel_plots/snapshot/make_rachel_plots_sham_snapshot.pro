dir = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla_v2/sham/'
outdir = dir+'plot_info_099/'
gfname = 'SHAM_catalog_099.fit_galaxies.fit'
hfname = 'SHAM_catalog_099.fit_halos.fit'
masscut = 1e12

make_rachel_plotdata_snapshot, outdir, dir, gfname=gfname, hfname=hfname, $
        masscut=masscut

end
