pro read_data_des_slice, num, hfile, format, $
	path=path, nfiles=nfiles, name = name, $
	skip_dr8_cat = skip_dr8_cat, read_full_files=read_full_files, $
	skip_stripe82_cat=skip_stripe82_cat, skip_photoz = skip_photoz, $
	skip_other_mags = skip_other_mags, first = first, $
	zmax = zmax, zmin=zmin, $
	omegam=omegam, omegal=omegal,$
	skip_all_galaxies=skip_all_galaxiees, skip_des_cat=skip_des_cat,$
	sva1=sva1, boxsize=boxsize

IF NOT KEYWORD_SET(omegam) then omegam = 0.25
IF NOT KEYWORD_SET(omegal) then omegal = 1.0 - omegam
IF NOT KEYWORD_SET(path) then path = './'
path += '/individual_slices/'
spawn, 'mkdir -p '+path
if not KEYWORD_SET(name) then name = 'PO_Carmen02_'
sigma_cut = 5.0
if not KEYWORD_SET(boxsize) then boxsize = 3000.
id_offset = 1000000
min_mass = 1e10

;;;read in all the data files
file = num+'/idl/'+name+'000_galaxies.fit'
g = mrdfits(file,1)
g.id = lindgen(N_ELEMENTS(g))+long64(num)*id_offset
file = num+'/idl/'+name+'000_sdss25.fit'
sdss = mrdfits(file,1)
file = num+'/idl/'+name+'000_deep.fit'
deep = mrdfits(file,1)
file = num+'/idl/'+name+'000_des.fit'
des = mrdfits(file,1)
file = num+'/idl/'+name+'000_johnson.fit'
johnson = mrdfits(file,1)
file = num+'/idl/'+name+'000_vista.fit'
vista = mrdfits(file,1)
file = num+'/idl/'+name+'000_flamex.fit'
flamex = mrdfits(file,1)
file = num+'/idl/'+name+'000_cfhtls.fit'
cfhtls = mrdfits(file,1)
file = num+'/idl/'+name+'000_euclid.fit'
euclid = mrdfits(file,1)
file = num+'/idl/'+name+'000_irac.fit'
irac = mrdfits(file,1)
file = num+'/idl/'+name+'000_wise.fit'
wise = mrdfits(file,1)
file = num+'/idl/'+name+'000_hsc.fit'
hsc = mrdfits(file,1)
file = num+'/idl/'+name+'000_lsst.fit'
lsst = mrdfits(file,1)
file = num+'/idl/'+name+'000_wfirst.fit'
wfirst = mrdfits(file,1)
  
;;;cut duplicate galaxies at halo centers
make_galaxies_unique, g, ii
g = g[ii]
sdss = sdss[ii]
deep = deep[ii]
des = des[ii]
vista = vista[ii]
johnson = johnson[ii]
flamex = flamex[ii]
cfhtls = cfhtls[ii]

;;option cut on redshift
if (keyword_set(zmin) and keyword_set(zmax)) then begin
  ii = where(g.z le zmax and g.z gt zmin)
  g = g[ii]
  sdss = sdss[ii]
  deep = deep[ii]
  des = des[ii]
  vista = vista[ii]
  johnson = johnson[ii]
  flamex = flamex[ii]
  cfhtls = cfhtls[ii]
endif  


;;;read in our halo file
case strupcase(format) of
  'FITS':begin
    hin = mrdfits(hfile,1)
  end
  'HLIST':begin
    read_hlist_general, hfile, hin
    cc = where(hin.pid lt 0)
    hin[cc].pid = hin[cc].id
  end
  'PARENTS':begin
    read_rockstar_parents, hfile, hin
    cc = where(hin.pid lt 0)
    hin[cc].pid = hin[cc].id
  end
endcase
id = where(hin.mvir ge min_mass)
hin = hin[ii]
h = fill_halo_tags(hin, omegam=omegam, omegal=omegal)

ng = N_ELEMENTS(g)
g.id = lindgen(ng)

;;;add all of our photometric errors
if not KEYWORD_SET(skip_other_mags) then begin

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
  mwrfits, cfhtls, path+'Mock_CFHTLS_mag_'+num+'.fit', /create
  
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
  mwrfits, rcs, path+'Mock_RCS_mag_'+num+'.fit', /create
  cfhtls = 0
  rcs = 0
  
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
  mwrfits, bcs, path+'Mock_BCS_mag_'+num+'.fit', /create
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
  mwrfits, deep, path+'Mock_DEEP2_mag_'+num+'.fit', /create
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
  mwrfits, vhs, path+'Mock_VHS_mag_'+num+'.fit', /create
  
  print, "adding viking photometric errors...."
  mock_error_apply,'VIKING',vista.omag,flux,ivar,omag,omagerr
  add_tags, vista, ['id', 'tmag', 'lmag', 'omagerr', 'flux', 'ivar'], ['0L', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], viking
  viking.id = g.id
  viking.tmag = vista.omag
  viking.omag = omag
  viking.omagerr = omagerr
  viking.flux = flux
  viking.ivar = 1./ivar^2
  mwrfits, viking, path+'Mock_VIKING_mag_'+num+'.fit', /create
  viking = 0
  vista = 0
  
  
  print, "Adding flamex errors...."
  mock_error_apply,'NDWFS',flamex.omag[0:2],flux1,ivar1,omag1,omagerr1
  mock_error_apply,'FLAMEX',flamex.omag[[3,5]],flux2,ivar2,omag2,omagerr2
  mock_error_apply,'IRAC',flamex.omag[6:7],flux3,ivar3,omag3,omagerr3
  
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
  mwrfits, ndwfs, path+'Mock_NDWFS_mag_'+num+'.fit', /create
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
  mwrfits, flamex_out, path+'Mock_FLAMEX_mag_'+num+'.fit', /create
  flamex_out = 0
  
  irac1 = create_struct('id', 0L, 'tmag', fltarr(2), 'lmag', fltarr(2), 'amag', fltarr(2), $
                        'omag', fltarr(2), 'omagerr', fltarr(2), $
                        'flux',fltarr(2), 'ivar', fltarr(2))
  irac = replicate(irac1,ng)
  irac.id = g.id
  irac.amag = flamex.amag[6:7]
  irac.tmag = flamex.omag[6:7]
  ;irac.lmag = flamex.lmag[6:7]
  irac.omag = omag3
  irac.omagerr= omagerr3
  irac.flux = flux3
  irac.ivar = 1./ivar3^2
  mwrfits, irac, path+'Mock_IRAC_mag_'+num+'.fit', /create
  irac = 0
  flamex = 0

  print "Creating WISE bands without errors..."
  wise1 = reate_struct('id', 0L, 'tmag', fltarr(4), 'lmag', fltarr(4), 'amag', fltarr(4), $
                          'omag', fltarr(4), 'omagerr', fltarr(4), $
                          'flux',fltarr(4), 'ivar', fltarr(4))
  wise_out = replicate(wise1,ng)
  wise_out.id= g.id
  wise_out.amag = wise.amag
  wise_out.omag = wise.omag
  if (tag_exist(wise, 'tmag')) then wise_out.tmag = wise.tmag
  if (tag_exist(wise, 'lmag')) then wise_out.lmag = wise.lmag
  mwrfits, wise_out, path+'Mock_WISE_mag_'+num+'.fit', /create
  wise_out = 0
  wise = 0

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
  mwrfits, dr8_out, path+'Mock_DR8_mag_'+num+'.fit', /create
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
  h_dr8 = h
  add_to_hod_nonperiodic, h_dr8, dr8, mass_def = 'mvir', omegam=omegam, omegal=omegal
  mwrfits, dr8, path+'DR8_Mock_obs_'+num+'.fit', /create
  mwrfits, g[ii], path+'DR8_Mock_galaxies_'+num+'.fit', /create
  mwrfits, sdss[ii], path+'DR8_Mock_sdss25_'+num+'.fit', /create
  mwrfits, h_dr8, path+'DR8_Mock_halos_'+num+'.fit', /create
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
mwrfits, s82_out, path+'Mock_Stripe82_mag_'+num+'.fit', /create
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
  h_stripe82 = h
  add_to_hod_nonperiodic, h_stripe82, stripe82_out, mass_def = 'mvir', omegam=omegam, omegal=omegal
  mwrfits, stripe82_out, path+'stripe82_Mock_obs_'+num+'.fit', /create
  mwrfits, g[ii], path+'stripe82_Mock_galaxies_'+num+'.fit', /create
  mwrfits, sdss[ii], path+'stripe82_Mock_sdss25_'+num+'.fit', /create
  mwrfits, h_stripe82, path+'stripe82_Mock_halos_'+num+'.fit', /create
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
  mwrfits, des_mag, path+'Mock_DES_mag_'+num+'.fit', /create
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
  h_24 = h
  fix_hod_non_periodic, h_24, des_5year, mass_def = 'mvir', min_mass=1e12, omegam=omegam, omegal=omegal
  mwrfits, des_5year, path+'DES_Mock_5year_'+num+'.fit', /create
  mwrfits, g[ii], path+'DES_Mock_galaxies_'+num+'.fit', /create
  mwrfits, sdss[ii], path+'DES_Mock_sdss25_'+num+'.fit', /create
  mwrfits, h_24, path+'DES_Mock_halos_'+num+'.fit', /create
endif

if KEYWORD_SET(SVA1) then begin
  print, 'Doing SVA1 mock...'
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
    h_sva1 = h
    add_to_hod_nonperiodic, h_sva1, sva1, $
	mass_def = 'mvir', omegam=omegam, omegal=omegal
    mwrfits, sva1, path+'SVA1_Mock_'+num+'.fit', /create
    mwrfits, h_sva1, path+'SVA1_Mock_halos_'+num+'.fit', /create
  endif
endif

if not KEYWORD_SET(skip_all_galaxies) then begin
  print, 'Making final HOD catalog of all galaxies..."
  add_to_hod_nonperiodic, h, g, mass_def = 'mvir', omegam=omegam, omegal=omegal
  mwrfits, g, path+'Mock_galaxies_'+num+'.fit', /create
  mwrfits, h, path+'Mock_halos_'+num+'.fit', /create
endif

end

