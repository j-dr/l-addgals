pro make_bcc_file

args = command_line_args(count=ct)
;args = '0'
;path = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Lb1050_v1.0'
;fbase = 'PO_Aardvark_1050'
;outpath = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Aardvark/Catalog_v1.0/individual_box_files/'
;hfile = '/lustre/ki/pfs/mbusha/projects/Aardvark/Lb1050/Quartant1/rockstar/halos_1050.fit'

;healpix_num = '0'
;path = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla-1/Lb1050_v2.0/'
;fbase = 'PO_Chinchilla-1_1050'
;outpath = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla-1/Catalog_v2.0//individual_box_files'
;hfile = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb1050/Octants_01/rockstar/halos_1050.fit'

healpix_num = args[0]
path = args[1]
fbase = args[2]
outpath = args[3]+'/'
hfile = args[4]
n_files = args[5]
help, healpix_num
spawn, 'mkdir -p '+outpath

print, "Making the catalog for healpix number "+healpix_num+'.'
make_bcc_healpix_catalog_final, path, fbase, outpath, healpix_num, hfile, n_files, zmin = 0.90, zmax = 2.0, /no_lensing, /galz, mag_lim_delta=2.0

end
