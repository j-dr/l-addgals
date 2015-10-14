pro make_bcc_file

arg = command_line_args(count=ct)
healpix_num = arg[0]
help, healpix_num
path = arg[1]
fbase = arg[2]
outpath = arg[3]
hfile = arg[4]
n_files = long(arg[5])
zmin = float(arg[6])
zmax = float(arg[7])

;healpix_num='0'
;PATH='./'
;FBASE='PO_Aardvark_4000'
;OUTPATH='/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.6a/individual_box_files/'
;HFILE='/lustre/ki/pfs/mbusha/projects/Aardvark/Lb4000/Quartant1/rockstar/halos_4000.fit'
;N_FILES=2
;ZMIN=0.9
;ZMAX=2.0


;path = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Lb1050_v0.6a'
;fbase = 'PO_Aardvark_1050'
;outpath = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.6a/individual_box_files/'
;hfile = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.6a/halos_1050.fit'
;n_files = 8

print, "Making the catalog for healpix number "+healpix_num+'.'
make_bcc_healpix_catalog_final, path, fbase, outpath, healpix_num, hfile, n_files, zmin = zmin, zmax = zmax, /no_lensing, mag_lim_delta=2.0

end
