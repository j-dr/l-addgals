pro make_big_bcc_catalog, inpath, outpath, hfile, mag_lim=mag_lim, no_lensing=no_lensing, boxsize=boxsize, first=first, last=last, llbins=llbins, zmax=zmax, zmin=zmin, zpadding=zpadding

if not KEYWORD_SET(first) then first = 0
if not KEYWORD_SET(last) then last = 800

for i = first, last do begin
  healpix_num = strcompress(string(i), /remove_all)
  gfiles = list_with_path('*galaxies.'+healpix_num+'.fit', inpath)
  if(N_ELEMENTS(gfiles) eq 1) then continue
  print, "Compiling Pixel ", healpix_num
  make_bcc_healpix_catalog,inpath, outpath, healpix_num, hfile, mag_lim=mag_lim, no_lensing=no_lensing, boxsize=boxsize, llbins=llbins, zmax=zmax, zmin=zmin, zpadding=zpadding
endfor

end
