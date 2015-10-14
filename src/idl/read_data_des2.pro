pro read_data_des2, deep=deep, des=des, vista=vista, johnson=johnson, $
	path=path, nfiles=nfiles, strlen=strlen, name = name, $
	skip_dr8_cat = skip_dr8_cat, read_full_files=read_full_files, $
	skip_stripe82_cat=skip_stripe82_cat, skip_photoz = skip_photoz, $
	skip_other_mags = skip_other_mags, first = first, zmax = zmax, $
	hfile = hfile, format=format, omegam=omegam, omegal=omegal,$
	skip_all_galaxies=skip_all_galaxiees, skip_des_cat=skip_des_cat,$
	sva1=sva1, boxsize=boxsize

IF NOT KEYWORD_SET(omegam) then omegam = 0.25
IF NOT KEYWORD_SET(omegal) then omegal = 1.0 - omegam
IF NOT KEYWORD_SET(first) then first = 0L
IF NOT KEYWORD_SET(path) then path = './'
path += '/'
spawn, 'mkdir -p '+path
if not KEYWORD_SET(name) then name = 'PO_Carmen02_'
if not KEYWORD_SET(nfiles) then nfiles = 32
if not KEYWORD_SET(strlen) then strlen = 3
if not KEYWORD_SET(zmax) then zmax = 1.33
mmin = 24.2
mmin = 26.
sigma_cut = 5.0
if not KEYWORD_SET(boxsize) then boxsize = 3000.

if not KEYWORD_SET(read_full_files) then begin

  num = '00'+strcompress(string(first), /remove_all)
  if (strlen eq 1) then num = '0'
  file = num+'/idl/'+name+'000_galaxies.fit'
  g = mrdfits(file,1)
  ii = where(g.omag(2) le mmin or g.amag(2) le -19)
  g = g[ii]
  file = num+'/idl/'+name+'000_sdss25.fit'
  sdss = mrdfits(file,1)
  sdss = sdss[ii]
  file = num+'/idl/'+name+'000_deep.fit'
  deep = mrdfits(file,1)
  deep = deep[ii]
  file = num+'/idl/'+name+'000_des.fit'
  des = mrdfits(file,1)
  des = des[ii]
  file = num+'/idl/'+name+'000_johnson.fit'
  johnson = mrdfits(file,1)
  johnson = johnson[ii]
  file = num+'/idl/'+name+'000_vista.fit'
  vista = mrdfits(file,1)
  vista = vista[ii]
  file = num+'/idl/'+name+'000_flamex.fit'
  flamex = mrdfits(file,1)
  flamex = flamex[ii]
  file = num+'/idl/'+name+'000_cfhtls.fit'
  cfhtls = mrdfits(file,1)
  cfhtls = cfhtls[ii]
  file = num+'/idl/'+name+'000_euclid.fit'
  euclid = mrdfits(file,1)
  euclid = euclid[ii]
  file = num+'/idl/'+name+'000_irac.fit'
  irac = mrdfits(file,1)
  irac = irac[ii]
  file = num+'/idl/'+name+'000_wise.fit'
  wise = mrdfits(file,1)
  wise = wise[ii]
  file = num+'/idl/'+name+'000_hsc.fit'
  hsc = mrdfits(file,1)
  hsc = hsc[ii]
  file = num+'/idl/'+name+'000_lsst.fit'
  lsst = mrdfits(file,1)
  lsst = lsst[ii]
  file = num+'/idl/'+name+'000_wfirst.fit'
  wfirst = mrdfits(file,1)
  wfirst = wfirst[ii]
  
  
  openw,5,'redshift_splits.txt'
  for i = first+1, nfiles - 1 do begin
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
     file = dir+'/idl/'+name+'000_euclid.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_euclid.fit'
     teuclid = mrdfits(file,1)
     teuclid = teuclid[ii]
     euclid = [euclid[ii0], teuclid[ii1]]
     teuclid = 0
     file = dir+'/idl/'+name+'000_irac.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_irac.fit'
     tirac = mrdfits(file,1)
     tirac = tirac[ii]
     irac = [irac[ii0], tirac[ii1]]
     tirac = 0
     file = dir+'/idl/'+name+'000_wise.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_wise.fit'
     twise = mrdfits(file,1)
     twise = twise[ii]
     wise = [wise[ii0], twise[ii1]]
     twise = 0
     file = dir+'/idl/'+name+'000_hsc.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_hsc.fit'
     thsc = mrdfits(file,1)
     thsc = thsc[ii]
     chsc = [hsc[ii0], thsc[ii1]]
     thsc = 0
     file = dir+'/idl/'+name+'000_lsst.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_lsst.fit'
     tlsst = mrdfits(file,1)
     tlsst = tlsst[ii]
     lsst = [lsst[ii0], tlsst[ii1]]
     tlsst = 0
     file = dir+'/idl/'+name+'000_wfirst.fit'
     if (strlen eq 1) then file = dir+'/idl/'+name+dir+'_wfirst.fit'
     twfirst = mrdfits(file,1)
     twfirst = twfirst[ii]
     cwfirst = [wfirst[ii0], twfirst[ii1]]
     twfirst = 0
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
  if (strlen ne 1) then mwrfits, euclid, path+'euclid_full.fit', /create
  if (strlen ne 1) then mwrfits, irac, path+'irac_full.fit', /create
  if (strlen ne 1) then mwrfits, wise, path+'wise_full.fit', /create
  if (strlen ne 1) then mwrfits, hsc, path+'hsc_full.fit', /create
  if (strlen ne 1) then mwrfits, lsst, path+'lsst_full.fit', /create
  if (strlen ne 1) then mwrfits, wfirst, path+'wfirst_full.fit', /create

endif else begin
  g = mrdfits(path+'galaxies_full.fit',1)
  sdss = mrdfits(path+'sdss_full.fit',1)
  des = mrdfits(path+'des_full.fit',1)
  if not KEYWORD_SET(skip_other_mags) then begin
    deep = mrdfits(path+'deep_full.fit',1)
    vista = mrdfits(path+'vista_full.fit',1)
;    johnson = mrdfits(path+'johnson_full.fit',1)
    flamex = mrdfits(path+'flamex_full.fit',1)
    if (strlen ne 1) then cfhtls = mrdfits(path+'cfhtls_full.fit',1)
  endif
endelse

dr8_maglim = [22.12, 22.60, 22.29, 21.85, 20.32]
stripe82_maglim = dr8_maglim + 2.0

;hfile = '~/ki04/projects/LasDamas/Simulations/Carmen/02/analysis/LC/dc5/output/groups/halos_000_1e12'
if not KEYWORD_SET(hfile) then begin
  hfile = '/lustre/ki/pfs/mbusha/projects/LasDamas/Simulations/Carmen/02/analysis/LC/high_z_wedge/fixed/64_files/rockstar/bak/out_0'
  format = 'rockstar'
endif
ng = N_ELEMENTS(g)
g.id = lindgen(ng)

;;;add all of our photometric errors
if not KEYWORD_SET(skip_other_mags) then begin

  if (strlen ne 1) then begin
  print, "making cfhtls phtometric errors..."
  mock_error_apply, 'CFHTLS', cfhtls.omag, flux, ivar, omag, omagerr
  add_tags, cfhtls, ['id', 'tmag', 'omagerr', 'flux', 'ivar'], ['0L', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], mout
  mout.id = g.id
  mout.tmag = mout.omag
  mout.omag = omag
  mout.omagerr = omagerr
  mout.flux = flux
  mout.ivar = 1./(ivar)^2
  cfhtls = mout
  mwrfits, cfhtls, path+'Mock_CFHTLS_mag.fit', /create
  endif
  
  if (strlen ne 1) then begin
  print, "making rcs phtometric errors..."
  mock_error_apply, 'RCS', cfhtls.omag[1:4], flux, ivar, omag, omagerr
  rcs1 = create_struct('ID', 0L, 'omag', fltarr(4), 'amag', fltarr(4), 'tmag', fltarr(4), 'omagerr', fltarr(4), 'flux', fltarr(4), 'ivar', fltarr(4))
  rcs = replicate(rcs1, ng)
  rcs.id = g.id
  rcs.tmag = cfhtls.tmag[1:4]
  rcs.amag = cfhtls.amag[1:4]
  rcs.omag = omag
  rcs.omagerr = omagerr
  rcs.flux = flux
  rcs.ivar = 1./(ivar)^2
  mwrfits, rcs, path+'Mock_RCS_mag.fit', /create
  cfhtls = 0
  rcs = 0
  endif
  
  print, "making BCS phtometric errors..."
  mock_error_apply, 'BCS', g.omag[1:4], flux, ivar, omag, omagerr
  bcs1 = create_struct('ID', 0L, 'omag', fltarr(4), 'amag', fltarr(4), 'tmag', fltarr(4), 'omagerr', fltarr(4), 'flux', fltarr(4), 'ivar', fltarr(4))
  bcs = replicate(bcs1, ng)
  bcs.id = g.id
  bcs.tmag = g.omag[1:4]
  bcs.amag = g.amag[1:4]
  bcs.omag = omag
  bcs.omagerr = omagerr
  bcs.flux = flux
  bcs.ivar = 1./(ivar)^2
  mwrfits, bcs, path+'Mock_BCS_mag.fit', /create
  bcs = 0
  
  
  print, "making deep2 phtometric errors..."
  mock_error_apply, 'DEEP2', deep.omag, flux, ivar, omag, omagerr
  add_tags, deep, ['ID', 'tmag', 'omagerr', 'flux', 'ivar'], ['0L', 'fltarr(3)', 'fltarr(3)', 'fltarr(3)', 'fltarr(3)'], mout
  mout.id = g.id
  mout.tmag = mout.omag
  mout.omag = omag
  mout.omagerr = omagerr
  mout.flux = flux
  mout.ivar = 1./(ivar)^2
  deep = mout
  mwrfits, deep, path+'Mock_DEEP2_mag.fit', /create
  deep = 0
  
  print, 'adding VHS errors....'
  mock_error_apply, 'VHS', vista.omag[2:4], flux, ivar, omag, omagerr
  vhs1 = create_struct('ID', 0L, 'tmag', fltarr(3),'lmag',fltarr(3),'amag',fltarr(3), $
                            'omag', fltarr(3), 'omagerr', fltarr(3), $
                            'flux',fltarr(3), 'ivar', fltarr(3))
  vhs = replicate(vhs1,ng)
  vhs.id = g.id
  vhs.tmag = vista.omag[2:4]
  ;vhs.lmag = vista.lmag[2:4]
  vhs.amag = vista.amag[2:4]
  vhs.omag = omag
  vhs.omagerr = omagerr
  vhs.flux = flux
  vhs.ivar = 1./ivar^2
  mwrfits, vhs, path+'Mock_VHS_mag.fit', /create
  
  print, "adding viking photometric errors...."
  mock_error_apply,'VIKING',vista.omag,flux,ivar,omag,omagerr
  add_tags, vista, ['id', 'tmag', 'lmag', 'omagerr', 'flux', 'ivar'], ['0L', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], viking
  viking.id = g.id
  viking.tmag = vista.omag
  viking.omag = omag
  viking.omagerr = omagerr
  viking.flux = flux
  viking.ivar = 1./ivar^2
  mwrfits, viking, path+'Mock_VIKING_mag.fit', /create
  viking = 0
  vista = 0
  
  
  print, "Adding flamex errors...."
  mock_error_apply,'NDWFS',flamex.omag[0:2],flux1,ivar1,omag1,omagerr1
  mock_error_apply,'FLAMEX',flamex.omag[[3,5]],flux2,ivar2,omag2,omagerr2
  mock_error_apply,'WISE',flamex.omag[6:7],flux3,ivar3,omag3,omagerr3
  
  ndwfs1 = create_struct('id', 0L, 'tmag', fltarr(3), 'lmag', fltarr(3), 'amag', fltarr(3), $
                         'omag', fltarr(3), 'omagerr', fltarr(3), $
                         'flux',fltarr(3), 'ivar', fltarr(3))
  ndwfs = replicate(ndwfs1,ng)
  ndwfs.id = g.id
  ndwfs.amag = flamex.amag[0:2]
  ndwfs.tmag = flamex.omag[0:2]
  ;ndwfs.lmag = flamex.lmag[0:2]
  ndwfs.omag = omag1
  ndwfs.omagerr = omagerr1
  ndwfs.flux = flux1
  ndwfs.ivar = 1./ivar1^2
  mwrfits, ndwfs, path+'Mock_NDWFS_mag.fit', /create
  ndwfs = 0
  
  flamex1 = create_struct('id', 0L, 'tmag', fltarr(3), 'lmag', fltarr(3), 'amag', fltarr(3), $
                          'omag', fltarr(3), 'omagerr', fltarr(3), $
                          'flux',fltarr(3), 'ivar', fltarr(3))
  flamex_out = replicate(flamex1,ng)
  flamex_out.id = g.id
  flamex_out.amag = flamex.amag[3:5]
  flamex_out.tmag = flamex.omag[3:5]
  ;flamex_out.lmag = flamex.lmag[3:5]
  flamex_out.omag[[3,5]] = omag2
  flamex_out.omagerr[[3,5]] = omagerr2
  flamex_out.flux[[3,5]] = flux2
  flamex_out.ivar[[3,5]] = 1./ivar2^2
  mwrfits, flamex_out, path+'Mock_FLAMEX_mag.fit', /create
  flamex_out = 0
  
  wise1 = create_struct('id', 0L, 'tmag', fltarr(2), 'lmag', fltarr(2), 'amag', fltarr(2), $
                        'omag', fltarr(2), 'omagerr', fltarr(2), $
                        'flux',fltarr(2), 'ivar', fltarr(2))
  wise = replicate(wise1,ng)
  wise.id = g.id
  wise.amag = flamex.amag[6:7]
  wise.tmag = flamex.omag[6:7]
  ;wise.lmag = flamex.lmag[6:7]
  wise.omag = omag3
  wise.omagerr= omagerr3
  wise.flux = flux3
  wise.ivar = 1./ivar3^2
  mwrfits, wise, path+'Mock_WISE_mag.fit', /create
  wise = 0
  flamex = 0

endif

if not KEYWORD_SET(skip_dr8_cat) then begin
  print, "making dr8 phtometric errors..."
  mock_error_apply,'DR8',g.omag,flux,ivar,omag,omagerr
  
  dr8_out1 = create_struct('id', 0L, 'tmag', fltarr(5), 'lmag', fltarr(5), 'amag', fltarr(5), $
                            'omag', fltarr(5), 'omagerr', fltarr(5), $
                            'flux',fltarr(5), 'ivar', fltarr(5))
  dr8_out = replicate(dr8_out1,ng)
  dr8_out.id = g.id
  dr8_out.amag = g.amag
  dr8_out.tmag = g.omag
  dr8_out.omag = omag
  dr8_out.omagerr = omagerr
  dr8_out.flux = flux
  dr8_out.ivar = 1./(ivar)^2
  mwrfits, dr8_out, path+'Mock_DR8_mag.fit', /create
  dr8_out = 0

  print, "Assembling the DR8 catalog...."
  add_tags, g, ['tmag', 'omagerr', 'flux', 'ivar'], ['fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], dr8
  dr8.tmag = dr8.omag
  dr8.omag = omag
  dr8.omagerr = omagerr
  dr8.flux = flux
  dr8.ivar = 1./(ivar)^2	
;  ii = where(dr8.omag[0] lt dr8_maglim[0] or $
;             dr8.omag[1] lt dr8_maglim[1] or $
;             dr8.omag[2] lt dr8_maglim[2] or $
;             dr8.omag[3] lt dr8_maglim[3] or $
;             dr8.omag[4] lt dr8_maglim[4])
  ii = where(dr8.flux[0]*sqrt(dr8.ivar[0]) ge sigma_cut or $
             dr8.flux[1]*sqrt(dr8.ivar[1]) ge sigma_cut or $
             dr8.flux[2]*sqrt(dr8.ivar[2]) ge sigma_cut or $
             dr8.flux[3]*sqrt(dr8.ivar[3]) ge sigma_cut or $
             dr8.flux[4]*sqrt(dr8.ivar[4]) ge sigma_cut)
  dr8 = dr8[ii]
  fix_hod_non_periodic, h_dr8, dr8, hfile=hfile, format=format, boxsize=boxsize, mass_def = 'mvir', min_mass=1e10, omegam=omegam, omegal=omegal
  itrain = get_dr8_training_set(dr8)
  mwrfits, dr8, path+'DR8_Mock_obs.fit', /create
  mwrfits, g[ii], path+'DR8_Mock_galaxies.fit', /create
  mwrfits, sdss[ii], path+'DR8_Mock_sdss25.fit', /create
  mwrfits, h_dr8, path+'DR8_Mock_halos.fit', /create
  mwrfits, dr8[itrain], path+'DR8_Mock_training_set.fit', /create
  if not KEYWORD_SET(skip_photoz) then $
    make_carlos_photoz, path+'DR8_Mock_training_set.fit', $
  	path+'DR8_Mock_obs.fit', maglim = 21.8, $
  	ofile = path+'DR8_Mock_zCarlos.fit', nne=100, $
  	zmin = 0.0, zmax = 1.1, grid = '35', res = '35'
  dr8 = 0
endif
 
print, "making stripe82 phtometric errors..."
mock_error_apply,'STRIPE82',g.omag,flux,ivar,omag,omagerr
add_tags, g, ['tmag', 'omagerr', 'flux', 'ivar'], ['fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], stripe82
stripe82.tmag = stripe82.omag
stripe82.omag = omag
stripe82.omagerr = omagerr
stripe82.flux = flux
stripe82.ivar = 1./(ivar)^2

s82_out1 = create_struct('id', 0L, 'tmag', fltarr(5), 'lmag', fltarr(5), 'amag', fltarr(5), $
                          'omag', fltarr(5), 'omagerr', fltarr(5), $
                          'flux',fltarr(5), 'ivar', fltarr(5))
s82_out = replicate(s82_out1,ng)
s82_out.id = g.id
s82_out.amag = g.amag
s82_out.tmag = g.omag
s82_out.omag = omag
s82_out.omagerr = omagerr
s82_out.flux = flux
s82_out.ivar = 1./(ivar)^2
mwrfits, s82_out, path+'Mock_Stripe82_mag.fit', /create
s82_out = 0

if not KEYWORD_SET(skip_stripe82_cat) then begin
  print, "Assembling the Stripe82 catalog...."
  add_tags, g, ['tmag', 'omagerr', 'flux', 'ivar'], ['fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], stripe82
  stripe82.tmag = stripe82.omag
  stripe82.omag = omag
  stripe82.omagerr = omagerr
  stripe82.flux = flux
  stripe82.ivar = 1./(ivar)^2
;  ii = where(stripe82.omag[0] lt stripe82_maglim[0] or $
;             stripe82.omag[1] lt stripe82_maglim[1] or $
;             stripe82.omag[2] lt stripe82_maglim[2] or $
;             stripe82.omag[3] lt stripe82_maglim[3] or $
;             stripe82.omag[4] lt stripe82_maglim[4])
  ii = where(stripe82.flux[0]*sqrt(stripe82.ivar[0]) ge sigma_cut or $
             stripe82.flux[1]*sqrt(stripe82.ivar[1]) ge sigma_cut or $
             stripe82.flux[2]*sqrt(stripe82.ivar[2]) ge sigma_cut or $
             stripe82.flux[3]*sqrt(stripe82.ivar[3]) ge sigma_cut or $
             stripe82.flux[4]*sqrt(stripe82.ivar[4]) ge sigma_cut)
  stripe82_out = stripe82[ii]
  fix_hod_non_periodic, h_stripe82, stripe82_out, hfile=hfile, format=format, boxsize=boxsize, mass_def = 'mvir', min_mass=1e12, omegam=omegam, omegal=omegal
  mwrfits, stripe82_out, path+'stripe82_Mock_obs.fit', /create
  mwrfits, g[ii], path+'stripe82_Mock_galaxies.fit', /create
  mwrfits, sdss[ii], path+'stripe82_Mock_sdss25.fit', /create
  mwrfits, h_stripe82, path+'stripe82_Mock_halos.fit', /create
  stripe82_out = 0
endif

if not KEYWORD_SET(skip_des_cat) then begin
  print, "making DES phtometric errors..."
  mock_error_apply,'DES',des.omag,flux,ivar,omag,omagerr
  des_mag1 = create_struct('id', 0L, 'tmag', fltarr(5), 'lmag', fltarr(5), 'amag', fltarr(5), $
                            'omag', fltarr(5), 'omagerr', fltarr(5), $
                            'flux',fltarr(5), 'ivar', fltarr(5))
  des_mag = replicate(des_mag1,ng)
  des_mag.id = g.id
  des_mag.amag = des.amag
  des_mag.tmag = des.omag
  des_mag.omag = omag
  des_mag.omagerr = omagerr
  des_mag.flux = flux
  des_mag.ivar = 1./(ivar)^2
  mwrfits, des_mag, path+'Mock_DES_mag.fit', /create
  des_mag = 0
  
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
  des_5year.haloid = -1
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
  des_5year.px = g.px
  des_5year.py = g.py
  des_5year.pz = g.pz
  des_5year.vx = g.vx
  des_5year.vy = g.vy
  des_5year.vz = g.vz
  des = 0
  
  des_maglim = [25.4, 24.9, 25.3, 24.6, 22.1] 
  ;ii = where(des_5year.omag[0] lt des_maglim[0] or $
  ;           des_5year.omag[1] lt des_maglim[1] or $
  ;           des_5year.omag[2] lt des_maglim[2] or $
  ;           des_5year.omag[3] lt des_maglim[3] or $
  ;           des_5year.omag[4] lt des_maglim[4])
    ii = where(des_5year.flux[0]*sqrt(des_5year.ivar[0]) ge sigma_cut or $
               des_5year.flux[1]*sqrt(des_5year.ivar[1]) ge sigma_cut or $
               des_5year.flux[2]*sqrt(des_5year.ivar[2]) ge sigma_cut or $
               des_5year.flux[3]*sqrt(des_5year.ivar[3]) ge sigma_cut or $
               des_5year.flux[4]*sqrt(des_5year.ivar[4]) ge sigma_cut)
  des_5year = des_5year[ii]
  fix_hod_non_periodic, h_24, des_5year, hfile=hfile, format=format, boxsize=boxsize, mass_def = 'mvir', min_mass=1e12, omegam=omegam, omegal=omegal
  hh = where(h_24.z le zmax or h_24.ngals gt 0)
  h_24 = h_24[hh]
  mwrfits, des_5year, path+'DES_Mock_5year.fit', /create
  mwrfits, g[ii], path+'DES_Mock_galaxies.fit', /create
  mwrfits, sdss[ii], path+'DES_Mock_sdss25.fit', /create
  mwrfits, h_24, path+'DES_Mock_halos.fit', /create
endif

if KEYWORD_SET(SVA1) then begin
  ;;;read in the mask/limmag files
  pixres = 9
  limmag_file = '/nfs/slac/g/ki/ki19/des/erykoff/des/sva1/downloads/sva1_coadd_limmags_pix9_dm.fit'
  limmag = mrdfits(limmag_file,1)
  maskfile = '/nfs/slac/g/ki/ki19/des/erykoff/des/sva1/downloads/sva1_coadd_maskstr_200.fit'
  mask = mrdfits(maskfile,1)

  ;;;cut the catalog based on the mask
  radec_to_simplepix,g.ra,g.dec,mask.pixres,pix
  calclambda_match_multi,mask.pix,pix,subb
  if (subb[0] gt 0) then begin
    sva1_g=g[subb]
    sva1_des = des[subb]
    sva1 = fill_bcc_truth(sva1_g, sva1_des)

    ;;;calcualte the limmag for each galaxy
    radec_to_simplepix, sva1.ra, sva1.dec, pixres, pix
    calclambda_match_multi, limmag.pixnum, pix, subb, suba=suba
    limmag = limmag[suba]
    sva1 = sva1[subb]

    ;;;apply the errors
    mock_error_apply_mask, sva1.tmag, limmag, omag, omagerr, flux, fluxerr
    sva1.omag = omag
    sva1.oamgerr = omagerr
    sva1.flux = flux
    sva1.ivar = 1./fluxerr^2

    ;;;cut the catalog based on errors
    ii = where(flux[0]/fluxerr[0] ge sva1_sigma_cut or $
	       flux[1]/fluxerr[1] ge sva1_sigma_cut or $
	       flux[2]/fluxerr[2] ge sva1_sigma_cut or $
	       flux[3]/fluxerr[3] ge sva1_sigma_cut or $
	       flux[4]/fluxerr[4] ge sva1_sigma_cut, count)
    sva1 = sva1[ii] 
    
    ;;;make the HOD information
    fix_hod_non_periodic, h_sva1, sva1, hfile=hfile, format=format, $
	boxsize=boxsize, mass_def = 'mvir', min_mass=1e12, omegam=omegam, $
	omegal=omegal
    hh = where(h_sva1.z le zmax or h_sva1.ngals gt 0)
    h_sva1 = h_sva1[hh]
    mwrfits, sva1, path+'SVA1_Mock.fit', /create
    mwrfits, h_sva1, path+'SVA1_Mock_halos.fit', /create
  endif
endif

if not KEYWORD_SET(skip_all_galaxies) then begin
  print, 'Making final HOD catalog of all galaxies..."
  fix_hod_non_periodic, h, g, hfile=hfile, format=format, boxsize=boxsize, mass_def = 'mvir', min_mass=1e12, omegam=omegam, omegal=omegal
  mwrfits, g, path+'Mock_galaxies.fit', /create
  mwrfits, h, path+'Mock_halos.fit', /create
endif

end

