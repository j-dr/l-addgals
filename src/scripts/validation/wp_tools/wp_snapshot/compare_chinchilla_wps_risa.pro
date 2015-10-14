plotfile = '~/temp/Chinchilla_wps.ps'
;labels = ['bcgs_1e12', 'bcgs_1e13', 'bcgs_1e14', 'no_bcgs']
labels = ['ADDGALS', 'SDSS']

sham_dir = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla_v2/sham/wp_info/'
model_dir = '/nfs/slac/g/ki/ki23/des/mbusha/Simulations/Chinchilla-tuning/downsample/addgals/single_epoch_fit/bcgs_1e12/idl/plot_info/'
;model_dir = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla-1/Catalog_v2.0/wp_cat/'
sdss_dir = '~/temp/data/'
dirs = [sham_dir, [model_dir], [sdss_dir]]

begplot, plotfile, /color, xsize = 7, ysize = 6
compare_wp, dirs, label = ['SHAM', [labels]]
endplot

end
