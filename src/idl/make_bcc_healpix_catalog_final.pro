pro make_bcc_healpix_catalog_final, path, fbase, outpath, healpix_num, hfile, n_files, mag_lim=mag_lim, boxsize=boxsize, llbins=llbins, zmax=zmax, zmin=zmin, zpadding=zpadding, mag_lim_delta=mag_lim_delta, no_lensing=no_lensing, omegam=omegam, omegal=omegal, galz=galz

  print, " "
  print, "  Making catalog for pixel ", healpix_num
  print, " "

  spawn, 'mkdir -p '+outpath

  if not KEYWORD_SET(zpadding) then zpadding = 0.04
  if not KEYWORD_SET(zmin) then zmin = 0.0
  if not KEYWORD_SET(zmax) then zmax = 100.0
  if not KEYWORD_SET(llbins) then llbins = 128
  id_offset = 10000000LL
  if not KEYWORD_SET(mag_lim_delta) then mag_lim_delta = 0.0

  ;;set our photometric depts
  if not KEYWORD_SET(mag_lim) then mag_lim = [25.5, 25.0, 24.4, 23.9, 22.0]
  
  print, 'zmin: ', zmin
  print, 'zmax: ', zmax

  base = path+'/'+healpix_num+'/000/idl/'+fbase+'_'+healpix_num+'.'
  tbase = base+'000_'
  print, 'tbase: ', tbase
  g = mrdfits(tbase+'galaxies.fit', 1)
  des = mrdfits(tbase+'des.fit',1)
  vista = mrdfits(tbase+'vista.fit',1)
  deep = mrdfits(tbase+'deep.fit',1)
  johnson = mrdfits(tbase+'johnson.fit',1)
  flamex = mrdfits(tbase+'flamex.fit',1)
  cfhtls = mrdfits(tbase+'cfhtls.fit',1)
  euclid = mrdfits(tbase+'euclid.fit',1)
  irac = mrdfits(tbase+'irac.fit',1)
  wise = mrdfits(tbase+'wise.fit',1)
  hsc = mrdfits(tbase+'hsc.fit',1)
  lsst = mrdfits(tbase+'lsst.fit',1)
  wfirst = mrdfits(tbase+'wfirst.fit',1)
  twomass = mrdfits(tbase+'twomass.fit',1)
  sdss25 = mrdfits(tbase+'sdss25.fit', 1)
  if KEYWORD_SET(galz) then  galz = mrdfits(tbase+'galz.fit', 1)
  ii = where(g.z ge zmin and g.z le zmax, count)
  print, 'size(ii):', size(ii)
  if (count eq 0) then ii = lindgen(N_ELEMENTS(g))
  g = g[ii]
  des = des[ii]
  vista = vista[ii]
  deep = deep[ii]
  johnson = johnson[ii]
  flamex = flamex[ii]
  cfhtls = cfhtls[ii]
  sdss25 = sdss25[ii]
  if KEYWORD_SET(galz) then galz = galz[ii]
  for i = 1, n_files - 1 do begin
    print, i
    num = '000'+strcompress(string(i), /remove_all)
    num = strmid(num, strlen(num)-3, 3)
    base = path+'/'+healpix_num+'/'+num+'/idl/'+fbase+'_'+healpix_num+'.'
    tbase = base+num+'_'
    print, 'tbase: ', tbase
    tg = mrdfits(tbase+'galaxies.fit', 1)
    tdes = mrdfits(tbase+'des.fit',1)
    tvista = mrdfits(tbase+'vista.fit',1)
    tdeep = mrdfits(tbase+'deep.fit',1)
    tjohnson = mrdfits(tbase+'johnson.fit',1)
    tflamex = mrdfits(tbase+'flamex.fit',1)
    tcfhtls = mrdfits(tbase+'cfhtls.fit',1)
    tsdss25 = mrdfits(tbase+'sdss25.fit', 1)
    teuclid = mrdfits(tbase+'euclid.fit',1)
    tirac = mrdfits(tbase+'irac.fit',1)
    twise = mrdfits(tbase+'wise.fit',1)
    thsc = mrdfits(tbase+'hsc.fit',1)
    tlsst = mrdfits(tbase+'lsst.fit',1)
    twfirst = mrdfits(tbase+'wfirst.fit',1)
    ttwomass = mrdfits(tbase+'twomass.fit',1)
    if KEYWORD_SET(galz) then tgalz = mrdfits(tbase+'galz.fit',1)
    tzmin = 0.5*(min(tg.z)+max(g.z))
    print, 'tzmin', tzmin
    print, 'min(tg.z):', min(tg.z)
    print, 'max(tg.z):', max(tg.z)
    print, 'max(g.z):', max(g.z)
    print, 'min(g.z):', min(g.z)
    ii1 = where(g.z le tzmin)
    print, 'min(ii1)', min(ii1)
    print, 'max(ii1)', max(ii1)
    print, 'size(g)', size(g)
    print, 'size(tg)', size(tg)
    ii2 = where(tg.z gt tzmin and tg.z le zmax, count)
    if (count eq 0) then continue
    g = [g[ii1],tg[ii2]]
    des = [des[ii1], tdes[ii2]]
    vista = [vista[ii1], tvista[ii2]]
    deep = [deep[ii1], tdeep[ii2]]
    johnson = [johnson[ii1], tjohnson[ii2]]
    flamex = [flamex[ii1], tflamex[ii2]]
    cfhtls = [cfhtls[ii1], tcfhtls[ii2]]
    sdss25 = [sdss25[ii1], tsdss25[ii2]]
    euclid = [euclid[ii1], teuclid[ii2]]
    irac = [irac[ii1], tirac[ii2]]
    wise = [wise[ii1], twise[ii2]]
    hsc = [hsc[ii1], thsc[ii2]]
    lsst = [lsst[ii1], tlsst[ii2]]
    wfirst = [wfirst[ii1], twfirst[ii2]]
    twomass = [twomass[ii1], ttwomass[ii2]]
    if KEYWORD_SET(galz) then galz = [galz[ii1], tgalz[ii2]]
  endfor
  
  print, 'done reading in files, now making redshift cuts'
  ;;;cut to a specified redshfit range
  ii = where(g.z ge zmin and g.z le zmax)
  g = g[ii]
  des = des[ii]
  vista = vista[ii]
  deep = deep[ii]
  johnson = johnson[ii]
  flamex = flamex[ii]
  cfhtls = cfhtls[ii]
  sdss25 = sdss25[ii]
  euclid = euclid[ii]
  irac = irac[ii]
  wise = wise[ii]
  hsc = hsc[ii]
  lsst = lsst[ii]
  wfirst = wfirst[ii]
  twomass = twomass[ii]
  if KEYWORD_SET(galz) then galz = galz[ii]
  print, 'done w/ redshift cuts, now making depth cuts'
  ;;;cut to DES depth, leaving some padding if specified
  ii = where(des.omag[0] le mag_lim[0] + mag_lim_delta OR $
	     des.omag[1] le mag_lim[1] + mag_lim_delta OR $
	     des.omag[2] le mag_lim[2] + mag_lim_delta OR $
	     des.omag[3] le mag_lim[3] + mag_lim_delta OR $
	     des.omag[4] le mag_lim[4] + mag_lim_delta)
  g = g[ii]
  des = des[ii]
  vista = vista[ii]
  deep = deep[ii]
  johnson = johnson[ii]
  flamex = flamex[ii]
  cfhtls = cfhtls[ii]
  sdss25 = sdss25[ii]
  euclid = euclid[ii]
  irac = irac[ii]
  wise = wise[ii]
  hsc = hsc[ii]
  lsst = lsst[ii]
  wfirst = wfirst[ii]
  twomass = twomass[ii]
  if KEYWORD_SET(galz) then galz = galz[ii]

  ng=n_elements(g)
  
  print, 'reading in cooper catalog'
  cooper_cat = mrdfits('/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/cooper/combined_dr6_cooper.fit',1)

  print, "Read the data, now saving alternate magnitude systems..."
  save_bcc_magnitudes, vista, outpath+fbase+'_vista_mag.'+healpix_num+'.fit', system='vista', /delete
  save_bcc_magnitudes, johnson, outpath+fbase+'_johnson_mag.'+healpix_num+'.fit', system='johnson', /delete
  save_bcc_magnitudes, deep, outpath+fbase+'_deep_mag.'+healpix_num+'.fit', system='deep', /delete
  save_bcc_magnitudes, flamex, outpath+fbase+'_flamex_mag.'+healpix_num+'.fit', system='flamex', /delete
  save_bcc_magnitudes, cfhtls, outpath+fbase+'_cfhtls_mag.'+healpix_num+'.fit', system='cfhtls', /delete
  save_bcc_magnitudes, euclid, outpath+fbase+'_euclid_mag.'+healpix_num+'.fit', system='euclid', /delete
  save_bcc_magnitudes, irac, outpath+fbase+'_irac_mag.'+healpix_num+'.fit', system='irac', /delete
  save_bcc_magnitudes, wise, outpath+fbase+'_wise_mag.'+healpix_num+'.fit', system='wise', /delete
  save_bcc_magnitudes, hsc, outpath+fbase+'_hsc_mag.'+healpix_num+'.fit', system='hsc', /delete
  save_bcc_magnitudes, lsst, outpath+fbase+'_lsst_mag.'+healpix_num+'.fit', system='lsst', /delete
  save_bcc_magnitudes, wfirst, outpath+fbase+'_wfirst_mag.'+healpix_num+'.fit', system='wfirst', /delete
  save_bcc_magnitudes, twomass, outpath+fbase+'_twomass_mag.'+healpix_num+'.fit', system='twomass', /delete
  if KEYWORD_SET(galz) then begin
     save_bcc_magnitudes, galz, outpath+fbase+'_galz.'+healpix_num+'.fit', system='galz', /delete
  end
  ;;;save the pure SDSS magnitudes (structure is different from above, and we don't want to delete them)
  sdss_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5))
  sdss_mag = replicate(sdss_mag1, ng)
  sdss_mag.amag = g.amag
  sdss_mag.tmag = g.omag
  sxaddpar, hdr_sdss, 'EXTNAME', 'SDSS'
  sxaddpar, hdr_sdss, 'IND_U',0
  sxaddpar, hdr_sdss, 'IND_G',1
  sxaddpar, hdr_sdss, 'IND_R',2
  sxaddpar, hdr_sdss, 'IND_I',3
  sxaddpar, hdr_sdss, 'IND_Z',4
  mwrfits, sdss_mag, outpath+fbase+'_sdss_mag.'+healpix_num+'.fit', hdr_sdss, /create

  ;;;make out main truth structure (called obs_full, for some archaic reason)
  obs_fullone = get_bcc_truth_structure()
  sdss_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5), 'omag', fltarr(5), 'omagerr', fltarr(5), 'flux', fltarr(5), 'ivar', fltarr(5))
  obs_full=replicate(obs_fullone,ng)
  sdss_mag = replicate(sdss_mag1, ng)

  
  ind = lindgen(ng)
  obs_full[ind].id=ind + long64(healpix_num)*id_offset
  obs_full[ind].ecatid = g.ecatid
  obs_full[ind].tra=g.ra
  obs_full[ind].tdec=g.dec
  obs_full[ind].z=g.z
  obs_full[ind].haloid=g.haloid
  obs_full[ind].rhalo=g.rhalo
  obs_full[ind].r200=g.r200
  obs_full[ind].m200=g.m200
  obs_full[ind].central=g.central
  obs_full[ind].px=g.px
  obs_full[ind].py=g.py
  obs_full[ind].pz=g.pz
  obs_full[ind].vx=g.vx
  obs_full[ind].vy=g.vy
  obs_full[ind].vz=g.vz
  obs_full[ind].mag_u = g.omag[0]
  obs_full[ind].tmag(0) = des.omag(0)
  obs_full[ind].tmag(1) = des.omag(1)
  obs_full[ind].tmag(2) = des.omag(2)
  obs_full[ind].tmag(3) = des.omag(3)
  obs_full[ind].tmag(4) = des.omag(4)
  obs_full[ind].amag(0) = des.amag(0)
  obs_full[ind].amag(1) = des.amag(1)
  obs_full[ind].amag(2) = des.amag(2)
  obs_full[ind].amag(3) = des.amag(3)
  obs_full[ind].amag(4) = des.amag(4)
  sdss_mag[ind].amag = g.amag
  sdss_mag[ind].tmag = g.omag
  obs_full[ind].ra=g.ra
  obs_full[ind].dec=g.dec
  obs_full[ind].epsilon = -99
  obs_full[ind].size = -99
  obs_full[ind].tsize = -99
  obs_full[ind].te = -99
  obs_full[ind].gamma1 = -99
  obs_full[ind].gamma2 = -99
  obs_full[ind].kappa = -99
  obs_full[ind].mu = -99
  obs_full[ind].lmag(0) = 99
  obs_full[ind].lmag(1) = 99
  obs_full[ind].lmag(2) = 99
  obs_full[ind].lmag(3) = 99
  obs_full[ind].lmag(4) = 99
  
  for i = 0, 4 do $
     obs_full.coeffs[i] = cooper_cat[obs_full.ecatid].coeffs[i]
     
  obs_full.index = lindgen(ng)
  print, "Finished filling truth catalogs."
  
  ;;;delete some variables to save memory
  cooper_cat = 0

  ;;;add photometric errors to the DES catalog
  mock_error_apply,'DES',obs_full.tmag,flux,ivar,omag,omagerr
  obs_full.omag = omag
  obs_full.omagerr = omagerr
  obs_full.flux = flux
  obs_full.ivar = 1./ivar^2

  ;;;add the HOD information and save the incrimental halo catalog
  print, "Adding HOD info...."
  h = mrdfits(hfile,1)
  print, "halofile read"
  add_to_hod_nonperiodic, h, obs_full, mag_ind = 1, mass_def = 'mvir', omegam=omegam, omegal=omegal
  print, "add_to_hod_nonperiodic done"
  g.haloid = obs_full.haloid
  g.ngals = obs_full.ngals
  g.m200 = obs_full.m200
  g.rhalo = obs_full.rhalo
  g.r200 = obs_full.r200
  hout = outpath+fbase+'_halos.'+healpix_num+'.fit'
  mwrfits, h, hout, /create
  print, "Wrote halos"
  print, 'making intermediate DR8 and Stripe82 errors.'
  ;;;add and save our intermediate DR8 photometry (for testing on the intermediate catalotgs)
  mock_error_apply, 'DR8', sdss_mag.tmag, flux, ivar, omag, omagerr
  sdss_mag.omag = omag
  sdss_mag.omagerr = omagerr
  sdss_mag.flux = flux
  sdss_mag.ivar = 1./ivar^2
  mwrfits, sdss_mag, outpath+fbase+'_DR8_Mock_mag.'+healpix_num+'.fit', /create

  ;;;add and save our intermediate DR8 photometry (for testing on the intermediate catalotgs)
  mock_error_apply, 'Stripe82', sdss_mag.tmag, flux, ivar, omag, omagerr
  sdss_mag.omag = omag
  sdss_mag.omagerr = omagerr
  sdss_mag.flux = flux
  sdss_mag.ivar = 1./ivar^2
  mwrfits, sdss_mag, outpath+fbase+'_Stripe82_Mock_mag.'+healpix_num+'.fit', /create

  print, "Finished adding photometric errors."

  for i = 0L, 4 do begin
     good = where(obs_full.omag[i] le 99, comp=bad)
     if (bad[0] ge 0) then obs_full[bad].omag[i] = 99
  endfor
    
  print, "Making visible catalog...."
  obs_visibleone = get_bcc_obs_structure()
  obs_visible = replicate(obs_visibleone, ng)

  print, "Filling the visible catalog..."
  obs_visible.id = obs_full.id
  obs_visible.index = obs_full.index
  obs_visible.ra = obs_full.ra
  obs_visible.dec = obs_full.dec
  obs_visible.mag_u = obs_full.mag_u
  obs_visible.mag_g = obs_full.omag[0]
  obs_visible.mag_r = obs_full.omag[1]
  obs_visible.mag_i = obs_full.omag[2]
  obs_visible.mag_z = obs_full.omag[3]
  obs_visible.mag_y = obs_full.omag[4]
  obs_visible.magerr_g = obs_full.omagerr[0]
  obs_visible.magerr_r = obs_full.omagerr[1]
  obs_visible.magerr_i = obs_full.omagerr[2]
  obs_visible.magerr_z = obs_full.omagerr[3]
  obs_visible.magerr_y = obs_full.omagerr[4]
  obs_visible.flux_g = obs_full.flux[0]
  obs_visible.flux_r = obs_full.flux[1]
  obs_visible.flux_i = obs_full.flux[2]
  obs_visible.flux_z = obs_full.flux[3]
  obs_visible.flux_y = obs_full.flux[4]
  obs_visible.ivar_g = obs_full.ivar[0]
  obs_visible.ivar_g = obs_full.ivar[0]
  obs_visible.ivar_g = obs_full.ivar[0]
  obs_visible.ivar_g = obs_full.ivar[0]
  obs_visible.epsilon1 = obs_full.epsilon[0]
  obs_visible.epsilon2 = obs_full.epsilon[1]
  obs_visible.size = obs_full.size
  zerr=0.03*randomn(seed,N_ELEMENTS(obs_full))
  obs_full.photoz_Gaussian=obs_full.z+zerr*(1.0+obs_full.z)
  obs_visible.photoz_Gaussian = obs_full.photoz_Gaussian
  
    ;;;define header information
    sxaddpar, hdr_t, 'EXTNAME', 'DES'
    sxaddpar, hdr_t, 'IND_G',0
    sxaddpar, hdr_t, 'IND_R',1
    sxaddpar, hdr_t, 'IND_I',2
    sxaddpar, hdr_t, 'IND_Z',3
    sxaddpar, hdr_t, 'IND_Y',4
    sxaddpar, hdr_t, 'LIMMAG_G', mag_lim[0]
    sxaddpar, hdr_t, 'LIMMAG_R', mag_lim[1]
    sxaddpar, hdr_t, 'LIMMAG_I', mag_lim[2]
    sxaddpar, hdr_t, 'LIMMAG_Z', mag_lim[3]
    sxaddpar, hdr_t, 'LIMMAG_Y', mag_lim[4]
    sxaddpar, hdr_t, 'AREA', 214.86
    sxaddpar, hdr_t, 'RA_MIN', 0.0
    sxaddpar, hdr_t, 'RA_MAX', 90.0
    sxaddpar, hdr_t, 'DEC_MIN', 0.0
    sxaddpar, hdr_t, 'DEC_MAX', 90.0
    sxaddpar, hdr_t, 'Z_LO', zmin
    sxaddpar, hdr_t, 'Z_HI', zmax
    sxaddpar, hdr_t, 'HEALPIX_NUM', long(healpix_num)

    sxaddpar, hdr_v, 'EXTNAME', 'DES'
    sxaddpar, hdr_v, 'IND_G',0
    sxaddpar, hdr_v, 'IND_R',1
    sxaddpar, hdr_v, 'IND_I',2
    sxaddpar, hdr_v, 'IND_Z',3
    sxaddpar, hdr_v, 'IND_Y',4
    sxaddpar, hdr_v, 'LIMMAG_G', mag_lim[0]
    sxaddpar, hdr_v, 'LIMMAG_R', mag_lim[1]
    sxaddpar, hdr_v, 'LIMMAG_I', mag_lim[2]
    sxaddpar, hdr_v, 'LIMMAG_Z', mag_lim[3]
    sxaddpar, hdr_v, 'LIMMAG_Y', mag_lim[4]
    sxaddpar, hdr_v, 'AREA', 220.552
    sxaddpar, hdr_v, 'RA_MIN', 10.0
    sxaddpar, hdr_v, 'RA_MAX', 30.0
    sxaddpar, hdr_v, 'DEC_MIN', 35.0
    sxaddpar, hdr_v, 'DEC_MAX', 50.0

  print, "Writing out first set of catalogs...."
;  mwrfits, obs_full[ii], outpath+fbase+'_truth_no_photoz.'+healpix_num+'.fit', hdr_t, /create
  mwrfits, obs_full, outpath+fbase+'_truth_no_photoz.'+healpix_num+'.fit', hdr_t, /create
  mwrfits, obs_visible, outpath+fbase+'_no_photoz.'+healpix_num+'.fit', hdr_v, /create

  print, "Finished combining the files."

end



