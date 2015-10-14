ps = 1

label = ['Fox_tuning_single_epoch']
plotfile = '~/temp/rdel_fox_'+label+'_099.ps'

dir = '/nfs/slac/g/ki/ki23/des/mbusha/Simulations/Chinchilla-tuning/downsample/addgals/single_epoch_fit/'+label+'/'
hv_path = dir+'hv_output/'
paramfile = dir+'NumericalParameters'
shamfile = '/nfs/slac/g/ki/ki21/cosmo/jderose/halos/sham_and_histograms/catalogs/Fox_v2/sham/SHAM_catalog_099.txt'


for i = 0L, N_ELEMENTS(label) - 1 do begin
  if (ps) then begplot, plotfile[i], /color, xsize = 7, ysize = 6

  check_rdel, hv_path=hv_path[i], paramfile=paramfile[i], shamfile=shamfile, simtype='SNAP', sham_rdel=sham_rdel, sham_mr=sham_mr

  if (ps) then endplot
endfor

end
