pro get_bossanova_galaxies, path, $
	truth_base=truth_base, sdss_base = sdss_base, $
	goutfile_a=goutfile_a, soutfile_a=soutfile_a, $
	goutfile_b=goutfile_b, soutfile_b=soutfile_b, $
	goutfile_c=goutfile_c, soutfile_c=soutfile_c, $
	goutfile_d=goutfile_d, soutfile_d=soutfile_d, $
	toutfile_a=toutfile_a, toutfile_b=toutfile_b, $
	toutfile_c=toutfile_c, toutfile_d=toutfile_d

if not KEYWORD_SET(truth_base) then truth_base = 'PO_Aardvark_truth'
if not KEYWORD_SET(sdss_base) then sdss_base = 'PO_Aardvark_sdss_mag'

old_training_file = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Aardvark/Catalog_v1.0/photoz_DR8_selected/Aardvark_v1.0_DR8_training_set_sdss_mag.fit'
t = mrdfits(old_training_file,1)

sq_deg = !PI
npatches_a = 15
npatches_b = 10
npatches_c = 10
npatches_d = 10

mlim_a = 20.75
mlim_b = 20.75
mlim_c = 21.0
mlim_d = 20.5
gmr_lim_a = 1.3
rmi_lim_a = 1.0
gmr_lim_b = 1.3
rmi_lim_b = 1.0
gmr_lim_c = 1.3
rmi_lim_c = 1.0
gmr_lim_d = 1.3
rmi_lim_d = 1.0


first = 208
nread = 0
for i = first, 400 do begin
  ;;;our files to read in
  file = path+'/truth/'+truth_base+'.'+strcompress(string(i),/remove_all)+'.fit'
  file2 = path+'/DR8/'+sdss_base+'.'+strcompress(string(i),/remove_all)+'.fit'
  if (file_test(file) eq 0) then continue
  g = mrdfits(file,1)
  s = mrdfits(file2, 1)
  sp = cut_to_dr8_training_limits(s, ii=ii)
  gp = g[ii]
  gmr = sp.omag[1] - sp.omag[2]
  rmi = sp.omag[2] - sp.omag[3]

  add_tags, sp, ['ra', 'dec', 'z'], ['0.', '0.', '0.'], sp2
  sp2.ra = gp.ra
  sp2.dec = gp.dec
  sp2.z = gp.z

  ;;;get the ra/dec range
  ramean = mean(gp.ra)
  decmean = mean(gp.dec)
  ddec = sqrt(sq_deg)
  dec1 = decmean - 0.5*ddec
  dec2 = dec1 + ddec
  theta1 = (90.-dec1)*!PI/180.
  theta2 = (90.-dec2)*!PI/180.
  dra = (sq_deg*(!PI/180.)^2/(cos(theta2) - cos(theta1)))*180./!PI
  ra1 = ramean - 0.5*dra
  ra2 = ramean + 0.5*dra

  ii_a = where(sp.omag[2] lt mlim_a and $
		gmr lt gmr_lim_a and $
		rmi lt rmi_lim_a and $
		gp.ra ge ra1 and gp.ra le ra2 and $
		gp.dec ge dec1 and gp.dec le dec2, tn_a)
  ii_b = where(sp.omag[2] lt mlim_b and $
                gmr lt gmr_lim_b and $
                rmi lt rmi_lim_b and $
                gp.ra ge ra1 and gp.ra le ra2 and $
                gp.dec ge dec1 and gp.dec le dec2, tn_b)
  ii_c = where(sp.omag[2] lt mlim_c and $
                gmr lt gmr_lim_c and $
                rmi lt rmi_lim_c and $
                gp.ra ge ra1 and gp.ra le ra2 and $
                gp.dec ge dec1 and gp.dec le dec2, tn_c)
  ii_d = where(sp.omag[2] lt mlim_d and $
                gmr lt gmr_lim_d and $
                rmi lt rmi_lim_d and $
                gp.ra ge ra1 and gp.ra le ra2 and $
                gp.dec ge dec1 and gp.dec le dec2, tn_d)
  print, tn_a, tn_b, tn_c, tn_d


  if (i eq first) then begin
    gspec_a = gp[ii_a]
    sspec_a = sp2[ii_a]
    if (nread lt npatches_b) then begin
      gspec_b = gp[ii_b]
      sspec_b = sp2[ii_b]
      gspec_c = gp[ii_c]
      sspec_c = sp2[ii_c]
      gspec_d = gp[ii_d]
      sspec_d = sp2[ii_d]
    endif
  endif else begin
    gspec_a = [gspec_a, gp[ii_a]]
    sspec_a = [sspec_a, sp2[ii_a]]
    if (nread lt npatches_b) then begin
      gspec_b = [gspec_b, gp[ii_b]]
      sspec_b = [sspec_b, sp2[ii_b]]
      gspec_c = [gspec_c, gp[ii_c]]
      sspec_c = [sspec_c, sp2[ii_c]]
      gspec_d = [gspec_d, gp[ii_d]]
      sspec_d = [sspec_d, sp2[ii_d]]
    endif
  endelse
  nread++
  if (nread ge npatches_a) then break
;asdf
endfor

mwrfits, gspec_a, goutfile_a, /create
mwrfits, sspec_a, soutfile_a, /create
mwrfits, [t, sspec_a], toutfile_a, /create
mwrfits, gspec_b, goutfile_b, /create
mwrfits, sspec_b, soutfile_b, /create
mwrfits, [t, sspec_b], toutfile_b, /create
mwrfits, gspec_c, goutfile_c, /create
mwrfits, sspec_c, soutfile_c, /create
mwrfits, [t, sspec_c], toutfile_c, /create
mwrfits, gspec_d, goutfile_d, /create
mwrfits, sspec_d, soutfile_d, /create
mwrfits, [t, sspec_d], toutfile_d, /create

end
