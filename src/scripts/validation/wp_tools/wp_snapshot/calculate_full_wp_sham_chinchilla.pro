num = '099'
rockfile = '/nfs/slac/g/ki/ki23/des/mbusha/Simulations/Chinchilla-tuning/analysis/hlists/hlist_'+num
shamfile = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla_v2/sham/SHAM_catalog_'+num+'.txt'
box_size = 400.
outpath = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla_v2/sham/wp_info/'

calculate_full_wp, box_size, outpath, rockfile=rockfile, shamfile=shamfile

end
