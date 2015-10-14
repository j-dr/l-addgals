pixelnum = lindgen(400)
pixelnum = strcompress(string(pixelnum),/remove_all)
npixels = N_ELEMENTS(pixelnum)

gfile_base = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Aardvark/Catalog_v1.0/photoz_tests/PRIMUS_XMM_Truth.'
sfile_base = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Aardvark/Catalog_v1.0/photoz_tests/PRIMUS_XMM_sdss_mag.'
outdir_base = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Aardvark/Catalog_v1.0/photoz_DR8_variance/'
f1 = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Aardvark/Catalog_v1.0/photoz_DR8/Aardvark_v1.0_DR8_training_set.fit'
f2 = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Aardvark/Catalog_v1.0/photoz_DR8/Aardvark_v1.0_DR8_training_set_sdss_mag.fit'

gp = mrdfits(f1,1)
sp = mrdfits(f2,1)

ii = where(gp.ra gt 157 and gp.ra lt 158.71    and gp.dec gt 4   and gp.dec lt 6, comp = jj)
g = gp[jj]
s = sp[jj]

for i = 0L, npixels - 1 do begin
  outdir = outdir_base+pixelnum[i]+'/photoz_files/'

  gfile = gfile_base+pixelnum[i]+'.fit'
  sfile = sfile_base+pixelnum[i]+'.fit'
  if (file_test(gfile) eq 0) then continue
  print, outdir
  spawn, "mkdir -p "+outdir

  tg = mrdfits(gfile,1)
  ts = mrdfits(sfile,1)

  gout = [g, tg]
  add_tags, ts, ['ra', 'dec', 'z'], ['0.', '0.', '0.'], ts2
  ts2.ra = tg.ra
  ts2.dec = tg.dec
  ts2.z = tg.z
  sout = [s, ts2]

  goutfile = outdir+'/Aardvark_v1.0_DR8_training_set_pixel_'+pixelnum[i]+'.fit'
  soutfile = outdir+'/Aardvark_v1.0_DR8_training_set_pixel_'+pixelnum[i]+'_sdss_mag.fit'
  mwrfits, gout, goutfile, /create
  mwrfits, sout, soutfile, /create
endfor

end
