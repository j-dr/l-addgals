plotfile = '~/temp/Fox_wps.ps'
labels = ['Fox_test_2']

sham_dir = '/nfs/slac/g/ki/ki21/cosmo/jderose/halos/sham_and_histograms/catalogs/Fox_v2/sham/wp_info/'
model_dir = '/nfs/slac/g/ki/ki22/cosmo/jderose/addgals/catalogs/'+labels+'/hv_output/wp_calc/'
dirs = [sham_dir, [model_dir]]

begplot, plotfile, /color, xsize = 7, ysize = 6
compare_wp, dirs, label = ['SHAM', [labels]]
endplot

end
