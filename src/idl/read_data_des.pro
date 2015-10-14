pro read_data_des, deep=deep, des=des, vista=vista, johnson=johnson, $
	path=path, nfiles=nfiles, strlen=strlen, name = name, $
	skip_dr8_cat = skip_dr8_cat, read_full_files=read_full_files, $
	skip_stripe82_cat=skip_stripe82_cat, skip_photoz = skip_photoz

IF NOT KEYWORD_SET(path) then path = './'
path += '/'
spawn, 'mkdir -p '+path
if not KEYWORD_SET(name) then name = 'PO_Carmen02_'
if not KEYWORD_SET(nfiles) then nfiles = 32
if not KEYWORD_SET(strlen) then strlen = 3
zmax = 1.33
mmin = 24.2
mmin = 26.

if not KEYWORD_SET(read_full_files) then begin

  num = '000'
  if (strlen eq 1) then num = '0'
  file = num+'/idl/'+name+num+'_galaxies.fit'
  g = mrdfits(file,1)
  ii = where(g.omag(2) le mmin or g.amag(2) le -19)
  g = g[ii]
  file = num+'/idl/'+name+num+'_sdss25.fit'
  sdss = mrdfits(file,1)
  sdss = sdss[ii]
  file = num+'/idl/'+name+num+'_deep.fit'
  deep = mrdfits(file,1)
  deep = deep[ii]
  file = num+'/idl/'+name+num+'_des.fit'
  des = mrdfits(file,1)
  des = des[ii]
  file = num+'/idl/'+name+num+'_johnson.fit'
  johnson = mrdfits(file,1)
  johnson = johnson[ii]
  file = num+'/idl/'+name+num+'_vista.fit'
  vista = mrdfits(file,1)
  vista = vista[ii]
  file = num+'/idl/'+name+num+'_flamex.fit'
  flamex = mrdfits(file,1)
  flamex = flamex[ii]
  file = num+'/idl/'+name+num+'_cfhtls.fit'
  cfhtls = mrdfits(file,1)
  cfhtls = cfhtls[ii]
  
  
  openw,5,'redshift_splits.txt'
  for i = 1, nfiles - 1 do begin
     print, ' '
     print, i
  ;   dir = strcompress(string(i), /remove_all)
     dir = '000'+strcompress(string(i), /remove_all)
     dir = strmid(dir, strlen(dir)-3,3)
     if (strlen eq 1) then dir = strcompress(string(i), /remove_all)
  ;   dir = bdir[i]
  ;   if (i eq 60) then dir = '../'+dir
     ;file = dir+'/idl/'+name+dir+'_galaxies.fit'  
     file = dir+'/idl/'+name+'000_galaxies.fit'  
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_galaxies.fit'  
     print, file
     print, ' '
     tg = mrdfits(file,1)
     ii = where(tg.omag(2) le mmin or tg.amag(2) lt -19)
     tg = tg[ii]
     range, tg.z
     range, g.z
     max0 = max(g.z)
     min0 = min(tg.z)
     edge = 0.5*(max0 + min0)
     printf, 5,edge
     ii0 = where(g.z le edge)
     ii1 = where(tg.z gt edge, n_new)
     if (n_new le 0) then continue
     g = [g[ii0],tg[ii1]]
     tg = 0
     file = dir+'/idl/'+name+'000_sdss25.fit'  
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_sdss25.fit'
     tsdss = mrdfits(file,1)
     tsdss = tsdss[ii]
     sdss = [sdss[ii0], tsdss[ii1]]
     file = dir+'/idl/'+name+'000_deep.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_deep.fit'
     tdeep = mrdfits(file,1)
     tdeep = tdeep[ii]
     deep = [deep[ii0], tdeep[ii1]]
     tdeep = 0
     file = dir+'/idl/'+name+'000_des.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_des.fit'
     tdes = mrdfits(file,1)
     tdes = tdes[ii]
     des = [des[ii0], tdes[ii1]]
     tdes = 0
     file = dir+'/idl/'+name+'000_vista.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_vista.fit'
     tvista = mrdfits(file,1)
     tvista = tvista[ii]
     vista = [vista[ii0], tvista[ii1]]
     tvista = 0
     file = dir+'/idl/'+name+'000_johnson.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_johnson.fit'
     tjohnson = mrdfits(file,1)
     tjohnson = tjohnson[ii]
     johnson = [johnson[ii0], tjohnson[ii1]]
     tjohnson = 0
     file = dir+'/idl/'+name+'000_flamex.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_flamex.fit'
     tflamex = mrdfits(file,1)
     tflamex = tflamex[ii]
     flamex = [flamex[ii0], tflamex[ii1]]
     tflamex = 0
     file = dir+'/idl/'+name+'000_cfhtls.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_cfhtls.fit'
     tcfhtls = mrdfits(file,1)
     tcfhtls = tcfhtls[ii]
     cfhtls = [cfhtls[ii0], tcfhtls[ii1]]
     tcfhtls = 0
  endfor
  close,5
  
  ii = where(g.z le zmax)
  g = g[ii]
  g.id = ii
  sdss = sdss[ii]
  deep = deep[ii]
  des = des[ii]
  vista = vista[ii]
  johnson = johnson[ii]
  flamex = flamex[ii]
  if (strlen ne 1) then cfhtls = cfhtls[ii]
  
  make_galaxies_unique, g, ii
  g = g[ii]
  g.id = ii
  sdss = sdss[ii]
  deep = deep[ii]
  des = des[ii]
  vista = vista[ii]
  johnson = johnson[ii]
  flamex = flamex[ii]
  if (strlen ne 1) then cfhtls = cfhtls[ii]
  
  mwrfits, g, path+'galaxies_full.fit', /create
  mwrfits, sdss, path+'sdss_full.fit', /create
  mwrfits, deep, path+'deep_full.fit', /create
  mwrfits, des, path+'des_full.fit', /create
  mwrfits, vista, path+'vista_full.fit', /create
  mwrfits, johnson, path+'johnson_full.fit', /create
  mwrfits, flamex, path+'flamex_full.fit', /create
  if (strlen ne 1) then mwrfits, cfhtls, path+'cfhtls_full.fit', /create

endif else begin
  g = mrdfits(path+'galaxies_full.fit',1)
  sdss = mrdfits(path+'sdss_full.fit',1)
  deep = mrdfits(path+'deep_full.fit',1)
  des = mrdfits(path+'des_full.fit',1)
  vista = mrdfits(path+'vista_full.fit',1)
  johnson = mrdfits(path+'johnson_full.fit',1)
  flamex = mrdfits(path+'flamex_full.fit',1)
  if (strlen ne 1) then cfhtls = mrdfits(path+'cfhtls_full.fit',1)
endelse

;hfile = '~/ki04/projects/LasDamas/Simulations/Carmen/02/analysis/LC/dc5/output/groups/halos_000_1e12'
hfile = '/lustre/ki/pfs/mbusha/projects/LasDamas/Simulations/Carmen/02/analysis/LC/high_z_wedge/fixed/64_files/rockstar/out_0'
ng = N_ELEMENTS(g)
g.id = lindgen(ng)

;;;add all of our photometric errors
if (strlen ne 1) then begin
print, "making cfhtls phtometric errors..."
mock_error_apply, 'CFHTLS', cfhtls.omag, flux, ivar, omag, omagerr
add_tags, cfhtls, ['tmag', 'omagerr', 'flux', 'ivar'], ['fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], mout
mout.tmag = mout.omag
mout.omag = omag
mout.omagerr = omagerr
mout.flux = flux
mout.ivar = 1./(ivar)^2
cfhtls = mout
endif

print, "making deep2 phtometric errors..."
mock_error_apply, 'DEEP2', deep.omag, flux, ivar, omag, omagerr
add_tags, deep, ['tmag', 'omagerr', 'flux', 'ivar'], ['fltarr(3)', 'fltarr(3)', 'fltarr(3)', 'fltarr(3)'], mout
mout.tmag = mout.omag
mout.omag = omag
mout.omagerr = omagerr
mout.flux = flux
mout.ivar = 1./(ivar)^2
deep = mout

if (strlen ne 1) then begin
print, "making rcs phtometric errors..."
mock_error_apply, 'RCS', cfhtls.omag[1:4], flux, ivar, omag, omagerr
rcs1 = create_struct('omag', fltarr(4), 'amag', fltarr(4), 'tmag', fltarr(4), 'omagerr', fltarr(4), 'flux', fltarr(4), 'ivar', fltarr(4))
rcs = replicate(rcs1, ng)
rcs.tmag = cfhtls.tmag[1:4]
rcs.amag = cfhtls.amag[1:4]
rcs.omag = omag
rcs.omagerr = omagerr
rcs.flux = flux
rcs.ivar = 1./(ivar)^2
endif

;add_dr8_photometric_errors_nanomaggies, g, dr8
print, "making dr8 phtometric errors..."
mock_error_apply,'DR8',g.omag,flux,ivar,omag,omagerr
add_tags, g, ['tmag', 'omagerr', 'flux', 'ivar'], ['fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], dr8
dr8.tmag = dr8.omag
dr8.omag = omag
dr8.omagerr = omagerr
dr8.flux = flux
dr8.ivar = 1./(ivar)^2

dr8_maglim = [22.12, 22.60, 22.29, 21.85, 20.32]
if not KEYWORD_SET(skip_dr8_cat) then begin
  ii = where(dr8.omag[0] lt dr8_maglim[0] or $
             dr8.omag[1] lt dr8_maglim[1] or $
             dr8.omag[2] lt dr8_maglim[2] or $
             dr8.omag[3] lt dr8_maglim[3] or $
             dr8.omag[4] lt dr8_maglim[4])
  ;fix_hod_non_periodic, h_dr8, dr8[ii], hfile=hfile, /ascii, boxsize=3000., llbins=128, /non_sham, /addgals
  fix_hod_non_periodic, h_dr8, dr8[ii], hfile=hfile, format='rockstar', boxsize=3000., mass_def = 'mvir'
  ;dr8.haloid = g.haloid
  ;dr8.rhalo = g.rhalo
  ;dr8.r200 = g.r200
  itrain = ii[get_dr8_training_set(dr8[ii])]
  mwrfits, dr8[ii], path+'DR8_Mock_obs.fit', /create
  mwrfits, g[ii], path+'DR8_Mock_galaxies.fit', /create
  mwrfits, des[ii], path+'DR8_Mock_des.fit', /create
  mwrfits, deep[ii], path+'DR8_Mock_deep.fit', /create
  mwrfits, vista[ii], path+'DR8_Mock_vista.fit', /create
  mwrfits, johnson[ii], path+'DR8_Mock_johnson.fit', /create
  mwrfits, flamex[ii], path+'DR8_Mock_flamex.fit', /create
  if (strlen ne 1) then mwrfits, cfhtls[ii], path+'DR8_Mock_cfhtls.fit', /create
  if (strlen ne 1) then mwrfits, rcs[ii], path+'DR8_Mock_rcs.fit', /create
  mwrfits, sdss[ii], path+'DR8_Mock_sdss25.fit', /create
  mwrfits, h_dr8, path+'DR8_Mock_halos.fit', /create
  mwrfits, dr8[itrain], path+'DR8_Mock_training_set.fit', /create
  if not KEYWORD_SET(skip_photoz) then $
    make_carlos_photoz, path+'DR8_Mock_training_set.fit', $
  	path+'DR8_Mock_obs.fit', maglim = 21.8, $
  	ofile = path+'DR8_Mock_zCarlos.fit', nne=100, $
  	zmin = 0.0, zmax = 1.1, grid = '35', res = '35'
endif
dr8_out1 = create_struct('omag', fltarr(5), 'amag', fltarr(5), 'tmag', fltarr(5), 'omagerr', fltarr(5), 'flux', fltarr(5), 'ivar', fltarr(5))
dr8_out = replicate(dr8_out1, ng)
dr8_out.omag = dr8.omag
dr8_out.omagerr = dr8.omagerr
dr8_out.flux = dr8.flux
dr8_out.ivar = dr8.ivar
dr8_out.tmag = dr8.tmag
dr8_out.amag = dr8.amag
dr8 = 0
 
print, "making stripe82 phtometric errors..."
mock_error_apply,'STRIPE82',g.omag,flux,ivar,omag,omagerr
add_tags, g, ['tmag', 'omagerr', 'flux', 'ivar'], ['fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], stripe82
stripe82.tmag = stripe82.omag
stripe82.omag = omag
omag = 0
stripe82.omagerr = omagerr
omagerr = 0
stripe82.flux = flux
flux = 0
stripe82.ivar = 1./(ivar)^2
ivar = 0
stripe82_maglim = dr8_maglim + 2.0
if not KEYWORD_SET(skip_stripe82_cat) then begin
  ii = where(stripe82.omag[0] lt stripe82_maglim[0] or $
             stripe82.omag[1] lt stripe82_maglim[1] or $
             stripe82.omag[2] lt stripe82_maglim[2] or $
             stripe82.omag[3] lt stripe82_maglim[3] or $
             stripe82.omag[4] lt stripe82_maglim[4])
  ;fix_hod_non_periodic, h_stripe82, stripe82[ii], hfile=hfile, /ascii, boxsize=3000., llbins=128, /non_sham, /addgals
  stripe82_out = stripe82[ii]
  fix_hod_non_periodic, h_stripe82, stripe82_out, hfile=hfile, format='rockstar', boxsize=3000., mass_def = 'mvir'
  mwrfits, stripe82_out, path+'stripe82_Mock_obs.fit', /create
  mwrfits, g[ii], path+'stripe82_Mock_galaxies.fit', /create
  mwrfits, des[ii], path+'stripe82_Mock_des.fit', /create
  mwrfits, deep[ii], path+'stripe82_Mock_deep.fit', /create
  mwrfits, vista[ii], path+'stripe82_Mock_vista.fit', /create
  mwrfits, johnson[ii], path+'stripe82_Mock_johnson.fit', /create
  mwrfits, flamex[ii], path+'stripe82_Mock_flamex.fit', /create
  if (strlen ne 1) then mwrfits, cfhtls[ii], path+'stripe82_Mock_cfhtls.fit', /create
  if (strlen ne 1) then mwrfits, rcs[ii], path+'stripe82_Mock_rcs.fit', /create
  mwrfits, sdss[ii], path+'stripe82_Mock_sdss25.fit', /create
  mwrfits, dr8_out[ii], path+'stripe82_Mock_dr8.fit', /create
  mwrfits, h_stripe82, path+'stripe82_Mock_halos.fit', /create
endif
s82_out = dr8_out
s82_out.omag = stripe82.omag
s82_out.omagerr = stripe82.omagerr
s82_out.flux = stripe82.flux
s82_out.ivar = stripe82.ivar
s82_out.tmag = stripe82.tmag
stripe82 = 0

print, "making DES phtometric errors..."
mock_error_apply,'DES',des.omag,flux,ivar,omag,omagerr
des_5year1 = get_bcc_truth_structure()
des_5year = replicate(des_5year1,ng)
des_5year.id = g.id
;des_5year.index = 
des_5year.ecatid = g.ecatid
if tag_exist(g, 'coeffs') then des_5year.coeffs = g.coeffs
des_5year.tmag = des.omag
des_5year.omag = omag
des_5year.flux = flux
des_5year.ivar = 1./ivar^2
des_5year.omagerr = omagerr
des_5year.amag = des.amag
des_5year.ra = g.ra
des_5year.dec = g.dec
des_5year.z = g.z
des_5year.haloid = g.haloid
des_5year.m200 = g.m200
des_5year.ngals = g.ngals
des_5year.r200 = g.r200
des_5year.central = g.central
des_5year.tra = g.ra
des_5year.tdec = g.dec
;des_5year.epsilon = 
;des_5year.gamma1 = 
;des_5year.gamma2 = 
;des_5year.kappa = 
;des_5year.mu = 
;des_5year.lmag = 
des_5year.mag_u = g.omag[0]
;des_5year.te = 
;des_5year.tsize = 
;des_5year.size = 
des_5year.photoz_Gaussian = g.z+0.03*randomn(seed,ng)*(1.0+g.z)
des_5year.px = g.vx
des_5year.py = g.vy
des_5year.pz = g.vz
des_5year.px = g.vx
des_5year.py = g.vy
des_5year.pz = g.vz


;add_des_photometric_errors_nanomaggies, g, des, des_5year
des_maglim = [25.4, 24.9, 25.3, 24.6, 22.1] 
ii = where(des_5year.omag[0] lt des_maglim[0] or $
           des_5year.omag[1] lt des_maglim[1] or $
           des_5year.omag[2] lt des_maglim[2] or $
           des_5year.omag[3] lt des_maglim[3] or $
           des_5year.omag[4] lt des_maglim[4])
;fix_hod_non_periodic, h_24, g[ii], hfile=hfile, /ascii, boxsize=3000., llbins=128, /non_sham, /addgals
fix_hod_non_periodic, h_24, des_5year[ii], hfile=hfile, format='rockstar', boxsize=3000., mass_def = 'mvir'
;des_5year[ii].haloid = g[ii].haloid
;des_5year[ii].rhalo = g[ii].rhalo
;des_5year[ii].r200 = g[ii].r200
hh = where(h_24.z le zmax or h_24.ngals gt 0)
h_24 = h_24[hh]
mwrfits, des_5year[ii], path+'DES_Mock_5year.fit', /create
mwrfits, g[ii], path+'DES_Mock_galaxies.fit', /create
mwrfits, des[ii], path+'DES_Mock_des.fit', /create
mwrfits, deep[ii], path+'DES_Mock_deep.fit', /create
mwrfits, vista[ii], path+'DES_Mock_vista.fit', /create
mwrfits, johnson[ii], path+'DES_Mock_johnson.fit', /create
mwrfits, flamex[ii], path+'DES_Mock_flamex.fit', /create
if (strlen ne 1) then mwrfits, cfhtls[ii], path+'DES_Mock_cfhtls.fit', /create
if (strlen ne 1) then mwrfits, rcs[ii], path+'DES_Mock_rcs.fit', /create
mwrfits, sdss[ii], path+'DES_Mock_sdss25.fit', /create
mwrfits, dr8[ii], path+'DES_Mock_DR8.fit', /create
mwrfits, s82_out[ii], path+'DES_Mock_Stripe82.fit', /create
mwrfits, h_24, path+'DES_Mock_halos.fit', /create

end

ii = get_training_set(g, random=0, nrandom = 3000, /stripe82, deep=deep, /vvds, /zcosmos)
srt = sort(g[ii].ra)
validation = ii[srt[0:19999]]
training = ii[srt[20000:N_ELEMENTS(ii)-1]]
mwrfits, des_5year[validation], path+'DES_Mock_validation_set.fit', /create
mwrfits, des_5year[training], path+'DES_Mock_training_set.fit', /create

;;make sub regions that are just 2 degrees wide
ramin = 10
dra = 2
for i = 0, 9 do begin
  tramin = ramin + dra*i
  tramax = tramin + dra
  ii = where(des_5year.ra ge tramin and des_5year.ra lt tramax)
  mwrfits, des_5year[ii], path+'DES_Mock_5year.'+strcompress(string(i),/remove_all)+'.fit', /create
  make_carlos_photoz, path+'/DES_Mock_training_set.fit', $
	path+'/DES_Mock_5year.'+strcompress(string(i),/remove_all)+'.fit', $
	ofile = path+'/DES_Mock_zCarlos.'+strcompress(string(i),/remove_all)
endfor


ii = where(g.omag[2] le 21)
fix_hod_non_periodic, h_21, g[ii], hfile=hfile, /ascii, boxsize=3000., llbins=128, /non_sham, /addgals
mwrfits, des_5year[ii], path+'DES_Mock_5year_21.fit', /create
mwrfits, g[ii], path+'DES_Mock_galaxies_21.fit', /create
mwrfits, des[ii], path+'DES_Mock_des_21.fit', /create
mwrfits, deep[ii], path+'DES_Mock_deep_21.fit', /create
mwrfits, vista[ii], path+'DES_Mock_vista_21.fit', /create
mwrfits, johnson[ii], path+'DES_Mock_johnson_21.fit', /create
mwrfits, flamex[ii], path+'DES_Mock_flamex_21.fit', /create
mwrfits, sdss[ii], path+'DES_Mock_sdss25_21.fit', /create
mwrfits, h_21, path+'DES_Mock_halos_21.fit', /create

;return, g
end
