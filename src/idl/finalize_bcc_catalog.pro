;pro finalize_bcc_catalog, inpath, outpath, shearpath, pixel, lensing
pro finalize_bcc_catalog

args=command_line_args(count=ct)
inpath = args[0]
outpath = args[1]
obase = args[2]
shearbase = args[3]
pixel = args[4]
lensing = long(args[5])
fbase = args[6]

nside_in = 4
nside_out = 8
sizes = ['1050', '2600', '4000']
;fbase = 'PO_Aardvark'
;inpath = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Aardvark/Catalog_v8c/'
;inpath = '/lustre/ki/pfs/mbusha/projects/Aardvark/catalogs/temp_lensing/'
;outpath = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Aardvark/Catalog_v8c/final_cat/'

combine_single_healpix_catalog, pixel, inpath, shearbase, fbase, outpath, sizes, nside_in, nside_out, /add_des_errors, /add_vista_errors, /add_dr8_errors, /add_stripe82_errors, /add_rcs_errors, add_lensing=lensing, obase=obase

end
