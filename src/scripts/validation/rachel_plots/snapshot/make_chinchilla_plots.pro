
dirlist = ['/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla_v2/sham/plot_info_099', $
;	'/nfs/slac/g/ki/ki23/des/mbusha/Simulations/Chinchilla-tuning/downsample/addgals/single_epoch_fit/no_bcgs/idl/plot_info/',$
;	'/nfs/slac/g/ki/ki23/des/mbusha/Simulations/Chinchilla-tuning/downsample/addgals/single_epoch_fit/bcgs_1e14/idl/plot_info/', $
	'/nfs/slac/g/ki/ki23/des/mbusha/Simulations/Chinchilla-tuning/downsample/addgals/single_epoch_fit/bcgs_1e13/idl/plot_info/', $
	'/nfs/slac/g/ki/ki23/des/mbusha/Simulations/Chinchilla-tuning/downsample/addgals/single_epoch_fit/bcgs_1e12/idl/plot_info/']

labels = ['SHAM', '1e13', '1e12']

outdir = '/nfs/slac/g/ki/ki23/des/mbusha/Simulations/Chinchilla-tuning/downsample/addgals/single_epoch_fit/rachel_plots/'

make_rachel_box_plots, dirlist, labels, outdir

end
