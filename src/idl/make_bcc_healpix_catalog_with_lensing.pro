pro make_bcc_healpix_catalog_final, path, fbase, outpath, healpix_num, hfile, n_files, mag_lim=mag_lim, no_lensing=no_lensing, boxsize=boxsize, llbins=llbins, zmax=zmax, zmin=zmin, zpadding=zpadding, mag_lim_delta=mag_lim_delta

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
  if not KEYWORD_SET(mag_lim) then mag_lim = [24.6, 24.1, 24.4, 23.8, 21.3]
  flux_sky = [4.6, 12.9, 17.7, 45.1, 14.9]          ;flux for sky
  zps = [26.4946, 26.5011, 26.2604, 25.9775, 23.7739] ;;;our photometric zero points
  apperture = 1.5
  pixel = 0.27
  apperture_area = !PI*(apperture/2.0)^2
  pixels_area = pixel^2
  npixels = apperture_area/pixels_area
  flux_gal = 10.^(0.4*(zps - mag_lim))
  exposure_time = (10./flux_gal)^2*(flux_sky*npixels + flux_gal)
  sig = 5
  flux_gal5 = (sig^2 + sqrt(sig^4 + 4*exposure_time*flux_sky*npixels*sig^2))/(2*exposure_time)
  mag_lim5 = zps - 2.5*alog10(flux_gal5)


  base = path+'/'+healpix_num+'/000/idl/'+fbase+'_'+healpix_num+'.'
  tbase = base+'000_'

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
  sdss25 = mrdfits(tbase+'sdss25.fit', 1)
  ii = where(g.z ge zmin and g.z le zmax, count)
  if (count eq 0) then ii = lindgen(N_ELEMENTS(g))
  g = g[ii]
  des = des[ii]
  vista = vista[ii]
  deep = deep[ii]
  johnson = johnson[ii]
  flamex = flamex[ii]
  cfhtls = cfhtls[ii]
  sdss25 = sdss25[ii]
  if NOT KEYWORD_SET(no_lensing) then begin
;    lg = mrdfits(path+'/from_joerg/galaxies_full_2.fit',1)
    lg = mrdfits(path+'/from_joerg/galaxies_full.fit',1)
    ldes = mrdfits(path+'/from_joerg/des_full.fit',1)
    lvista = mrdfits(path+'/from_joerg/vista_full.fit',1)
    ldeep = mrdfits(path+'/from_joerg/deep_full.fit',1)
    ljohnson = mrdfits(path+'/from_joerg/johnson_full.fit',1)
    lflamex = mrdfits(path+'/from_joerg/flamex_full.fit',1)
    lcfhtls = mrdfits(path+'/from_joerg/cfhtls_full.fit',1)
    leuclid = mrdfits(path+'/from_joerg/euclid_full.fit',1)
    lirac = mrdfits(path+'/from_joerg/irac_full.fit',1)
    lwise = mrdfits(path+'/from_joerg/wise_full.fit',1)
    lhsc = mrdfits(path+'/from_joerg/hsc_full.fit',1)
    llsst = mrdfits(path+'/from_joerg/lsst_full.fit',1)
    lwfirst = mrdfits(path+'/from_joerg/wfirst_full.fit',1)
  endif

  for i = 1, n_files - 1 do begin
    num = '000'+strcompress(string(i), /remove_all)
    num = strmid(num, strlen(num)-3, 3)
    base = path+'/'+healpix_num+'/'+num+'/idl/'+fbase+'_'+healpix_num+'.'
    tbase = base+num+'_'

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

    tzmin = 0.5*(min(tg.z)+max(g.z))
    ii1 = where(g.z le tzmin)
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
    irac = [irac[ii1], tirac[ii2]]
    wise = [wise[ii1], twise[ii2]]
    hsc = [hsc[ii1], thsc[ii2]]
    lsst = [lsst[ii1], tlsst[ii2]]
    wfirst = [wfirst[ii1], twfirst[ii2]]
  endfor

  ii = where(g.z ge zmin and g.z le zmax)
  g = g[ii]
  des = des[ii]
  vista = vista[ii]
  deep = deep[ii]
  johnson = johnson[ii]
  flamex = flamex[ii]
  cfhtls = cfhtls[ii]
  sdss25 = sdss25[ii]
  irac = irac[ii]
  wise = wise[ii]
  hsc = hsc[ii]
  lsst = lsst[ii]
  wfirst = wfirst[ii]

  ng=n_elements(g)
  if KEYWORD_SET(no_lensing) then ngnew = n_elements(g) else $
    ngnew = n_elements(lg)

  h = mrdfits(hfile,1)

  cooper_cat = mrdfits('/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/cooper/combined_dr6_cooper.fit',1)

  print, "Read the data, now making initial structures..."
  obs_fullone = get_bcc_truth_structure()
;  obs_fullone=create_struct('id',0LL,'index', 0L, 'ecatid', 0L, 'coeffs', fltarr(5), $
;                            'tmag',fltarr(5),'omag',fltarr(5), $
;			    'flux', fltarr(5), 'ivar', fltarr(5), $
;                            'omagerr',fltarr(5),'amag', fltarr(5), $
;                            'ra',0.0,'dec',0.0,'z',0.0,$
;                            'haloid',0L,'rhalo',0.0, 'm200', 0., 'ngals', 0L, $
;                            'r200',0.0,'central',0, 'tra', 0., 'tdec', 0., $
;                            'epsilon', fltarr(2), 'gamma1', 0., 'gamma2', 0., $
;                            'kappa', 0., 'mu', 0., 'lmag', fltarr(5), 'mag_u', 0.,$
;                            'te', fltarr(2), 'tsize', 0., 'size', 0., $
;                            'zCarlos',0.0, 'ArborZ', 0., $
;                            'ArborZ_err', 0., 'ANNZ', 0., 'ANNZ_err', 0.,$
;                            'photoz_Gaussian', 0., $
;                            'px', 0., 'py', 0., 'pz', 0., $
;                            'vx', 0., 'vy', 0., 'vz', 0.)
  sdss_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5), 'lmag', fltarr(5))
  vista_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5), 'lmag', fltarr(5))
  johnson_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5), 'lmag', fltarr(5))
  deep_mag1 = create_struct('amag', fltarr(3), 'tmag', fltarr(3), 'lmag', fltarr(3))
  flamex_mag1 = create_struct('amag', fltarr(8), 'tmag', fltarr(8), 'lmag', fltarr(8))
  cfhtls_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5), 'lmag', fltarr(5))
  euclid_mag1 = create_struct('amag', fltarr(4), 'tmag', fltarr(4), 'lmag', fltarr(4))
  irac_mag1 = create_struct('amag', fltarr(4), 'tmag', fltarr(4), 'lmag', fltarr(4))
  wise_mag1 = create_struct('amag', fltarr(4), 'tmag', fltarr(4), 'lmag', fltarr(4))
  hsc_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5), 'lmag', fltarr(5))
  lsst_mag1 = create_struct('amag', fltarr(8), 'tmag', fltarr(8), 'lmag', fltarr(8))
  wfirst_mag1 = create_struct('amag', fltarr(4), 'tmag', fltarr(4), 'lmag', fltarr(4))

  obs_full=replicate(obs_fullone,ngnew)
  sdss_mag = replicate(sdss_mag1, ngnew)
  vista_mag = replicate(vista_mag1, ngnew)
  johnson_mag = replicate(johnson_mag1, ngnew)
  deep_mag = replicate(deep_mag1, ngnew)
  flamex_mag = replicate(flamex_mag1, ngnew)
  cfhtls_mag = replicate(cfhtls_mag1, ngnew)
  euclid_mag = replicate(euclid_mag1, ngnew)
  irac_mag = replicate(irac_mag1, ngnew)
  wise_mag = replicate(wise_mag1, ngnew)
  hsc_mag = replicate(hsc_mag1, ngnew)
  lsst_mag = replicate(lsst_mag1, ngnew)
  wfirst_mag = replicate(wfirst_mag1, ngnew)

  ind = lindgen(ng)
  obs_full[ind].id=ind + float(healpix_num)*id_offset
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
  vista_mag[ind].amag = vista.amag
  vista_mag[ind].tmag = vista.omag
  johnson_mag[ind].amag = johnson.amag
  johnson_mag[ind].tmag = johnson.omag
  deep_mag[ind].amag = deep.amag
  deep_mag[ind].tmag = deep.omag
  flamex_mag[ind].amag = flamex.amag
  flamex_mag[ind].tmag = flamex.omag
  cfhtls_mag[ind].amag = cfhtls.amag
  cfhtls_mag[ind].tmag = cfhtls.omag
  euclid_mag[ind].amag = euclid.amag
  euclid_mag[ind].tmag = euclid.omag
  irac_mag[ind].amag = irac.amag
  irac_mag[ind].tmag = irac.omag
  wise_mag[ind].amag = wise.amag
  wise_mag[ind].tmag = wise.omag
  hsc_mag[ind].amag = hsc.amag
  hsc_mag[ind].tmag = hsc.omag
  lsst_mag[ind].amag = lsst.amag
  lsst_mag[ind].tmag = lsst.omag
  wfirst_mag[ind].amag = wfirst.amag
  wfirst_mag[ind].tmag = wfirst.omag


  if KEYWORD_SET(no_lensing) then begin
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
  endif else begin
    print, "Finding multiply-imaged galaxies...."
    for i = ng, ngnew - 1 do begin
      match = where(g.id eq lg[i].id, count)
      if (count ne 1) then print, 'Error with matching!!!!'
      obs_full[i].id=g[match].id
      obs_full[i].ecatid=g[match].ecatid
       obs_full[i].tra=g[match].ra
      obs_full[i].tdec=g[match].dec
      obs_full[i].z=g[match].z
      obs_full[i].haloid=g[match].haloid
      obs_full[i].rhalo=g[match].rhalo
      obs_full[i].r200=g[match].r200
      obs_full[i].m200=g[match].m200
      obs_full[i].central=g[match].central
      obs_full[i].px=g[match].px
      obs_full[i].py=g[match].py
      obs_full[i].pz=g[match].pz
      obs_full[i].vx=g[match].vx
      obs_full[i].vy=g[match].vy
      obs_full[i].vz=g[match].vz
      obs_full[i].mag_u = g[match].omag[0]
      obs_full[i].tmag(0) = des[match].omag(0)
      obs_full[i].tmag(1) = des[match].omag(1)
      obs_full[i].tmag(2) = des[match].omag(2)
      obs_full[i].tmag(3) = des[match].omag(3)
      obs_full[i].tmag(4) = des[match].omag(4)
      obs_full[i].amag(0) = des[match].amag(0)
      obs_full[i].amag(1) = des[match].amag(1)
      obs_full[i].amag(2) = des[match].amag(2)
      obs_full[i].amag(3) = des[match].amag(3)
      obs_full[i].amag(4) = des[match].amag(4)
      vista_mag[i].amag = vista[match].amag
      vista_mag[i].tmag = vista[match].omag
      johnson_mag[i].amag = johnson[match].amag
      johnson_mag[i].tmag = johnson[match].omag
      deep_mag[i].amag = deep[match].amag
      deep_mag[i].tmag = deep[match].omag
      flamex_mag[i].amag = flamex[match].amag
      flamex_mag[i].tmag = flamex[match].omag
      cfhtls_mag[i].amag = cfhtls[match].amag
      cfhtls_mag[i].tmag = cfhtls[match].omag
      sdss_mag[i].amag = g[match].amag
      sdss_mag[i].tmag = g[match].omag
      euclid_mag[i].amag = euclid[match].amag
      euclid_mag[i].tmag = euclid[match].omag
      irac_mag[i].amag = irac[match].amag
      irac_mag[i].tmag = irac[match].omag
      wise_mag[i].amag = wise[match].amag
      wise_mag[i].tmag = wise[match].omag
      hsc_mag[i].amag = hsc[match].amag
      hsc_mag[i].tmag = hsc[match].omag
      lsst_mag[i].amag = lsst[match].amag
      lsst_mag[i].tmag = lsst[match].omag
      wfirst_mag[i].amag = wfirst[match].amag
      wfirst_mag[i].tmag = wfirst[match].omag
      if KEYWORD_SET(no_lensing) then begin
        obs_full[i].ra=g[match].ra
        obs_full[i].dec=g[match].dec
        obs_full[i].epsilon = -99
        obs_full[i].size = -99
        obs_full[i].tsize = -99
        obs_full[i].te = -99
        obs_full[i].gamma1 = -99
        obs_full[i].gamma2 = -99
        obs_full[i].kappa = -99
        obs_full[i].mu = -99
        obs_full[i].lmag(0) = 99
        obs_full[i].lmag(1) = 99
        obs_full[i].lmag(2) = 99
        obs_full[i].lmag(3) = 99
        obs_full[i].lmag(4) = 99
      endif
    endfor
  endelse

  for i = 0, 4 do $
     obs_full.coeffs[i] = cooper_cat[obs_full.ecatid].coeffs[i]
     
  ng = ngnew
  obs_full.index = lindgen(ng)
  if NOT KEYWORD_SET(no_lensing) then begin
    obs_full.ra=lg.ora
    obs_full.dec=lg.odec
    obs_full.epsilon = lg.epsilon
    obs_full.size = lg.osize
    obs_full.tsize = lg.tsize
    obs_full.te = lg.te
    obs_full.gamma1 = lg.gamma1
    obs_full.gamma2 = lg.gamma2
    obs_full.kappa = lg.kappa
    obs_full.mu = lg.mu
    obs_full.lmag(0) = ldes.lmag(0)
    obs_full.lmag(1) = ldes.lmag(1)
    obs_full.lmag(2) = ldes.lmag(2)
    obs_full.lmag(3) = ldes.lmag(3)
    obs_full.lmag(4) = ldes.lmag(4)
    sdss_mag.lmag = lg.lmag
    vista_mag.lmag = lvista.lmag
    johnson_mag.lmag = ljohnson.lmag
    deep_mag.lmag = ldeep.lmag
    flamex_mag.lmag = lflamex.lmag
    cfhtls_mag.lmag = lcfhtls.lmag
    euclid_mag.lmag = euclid.lmag
    irac_mag.lmag = irac.lmag
    wise_mag.lmag = wise.lmag
    hsc_mag.lmag = hsc.lmag
    lsst_mag.lmag = lsst.lmag
    wfirst_mag.lmag = wfirst.lmag
  endif
  print, "Finished filling initial catalogs."
  
  ;;;delete some variables to save memory
  ldes = 0
  lg = 0
  lvista = 0
  ljohnson = 0
  ldeep = 0
  cooper_cat = 0

  if KEYWORD_SET(no_lensing) then add_des_photometric_errors,obs_full, ka, obs_full2, maglim=mag_lim, /no_new_struct
  if NOT KEYWORD_SET(no_lensing) then add_des_photometric_errors,obs_full, ka, obs_full2, maglim=mag_lim, /no_new_struct, /lmag
  obs_full = obs_full2

  print, "Adding HOD info...."
  add_to_hod_nonperiodic, h, obs_full, mag_ind = 1, mass_def = 'mvir'
  g.haloid = obs_full.haloid
  g.ngals = obs_full.ngals
  g.m200 = obs_full.m200
  g.rhalo = obs_full.rhalo
  g.r200 = obs_full.r200
  
;  hout = strmid(hfile, 0, strlen(hfile)-3) + healpix_num+'.fit'
  hout = outpath+fbase+'_halos.'+healpix_num+'.fit'
  mwrfits, h, hout, /create

  add_dr8_photometric_errors, g, dr8
  dr8_maglim = [22.12, 22.60, 22.29, 21.85, 20.32]
  ii = where(dr8.omag[0] lt dr8_maglim[0] or $
             dr8.omag[1] lt dr8_maglim[1] or $
             dr8.omag[2] lt dr8_maglim[2] or $
             dr8.omag[3] lt dr8_maglim[3] or $
             dr8.omag[4] lt dr8_maglim[4])
  mwrfits, dr8[ii], outpath+fbase+'_DR8_Mock_obs.'+healpix_num+'.fit', /create
  mwrfits, g[ii], outpath+fbase+'_DR8_Mock_galaxies.'+healpix_num+'.fit', /create

  add_stripe82_photometric_errors, g, stripe82
  stripe82_maglim = dr8_maglim + 2.0
  ii = where(stripe82.omag[0] lt stripe82_maglim[0] or $
             stripe82.omag[1] lt stripe82_maglim[1] or $
             stripe82.omag[2] lt stripe82_maglim[2] or $
             stripe82.omag[3] lt stripe82_maglim[3] or $
             stripe82.omag[4] lt stripe82_maglim[4])
  mwrfits, stripe82[ii], outpath+fbase+'_stripe82_Mock_obs.'+healpix_num+'.fit', /create
  mwrfits, g[ii], outpath+fbase+'_stripe82_Mock_galaxies.'+healpix_num+'.fit', /create



  print, "Finished adding photometric errors."

  print, "Cutting the catalog..."
;  ii = where(obs_full.omag[0] le mag_lim5[0]+mag_lim_delta or $
;             obs_full.omag[1] le mag_lim5[1]+mag_lim_delta or $
;             obs_full.omag[2] le mag_lim5[2]+mag_lim_delta or $
;             obs_full.omag[3] le mag_lim5[3]+mag_lim_delta or $
;             obs_full.omag[4] le mag_lim5[4]+mag_lim_delta, count)
  ii = where(obs_full.tmag[0] le mag_lim5[0]+mag_lim_delta or $
             obs_full.tmag[1] le mag_lim5[1]+mag_lim_delta or $
             obs_full.tmag[2] le mag_lim5[2]+mag_lim_delta or $
             obs_full.tmag[3] le mag_lim5[3]+mag_lim_delta or $
             obs_full.tmag[4] le mag_lim5[4]+mag_lim_delta, count)
  if (count gt 0) then begin
  ;obs_full = obs_full[ii]
  ;sdss_mag = sdss_mag[ii]
  ;vista_mag = vista_mag[ii]
  ;johnson_mag = johnson_mag[ii]
  ;deep_mag = deep_mag[ii]
  ;flamex_mag = flamex_mag[ii]
  ;cfhtls_mag = cfhtls_mag[ii]
    ng = count
    for i = 0L, 4 do begin
       good = where(obs_full.omag[i] le 99, comp=bad)
       if (bad[0] ge 0) then obs_full[bad].omag[i] = 99
    endfor
    
    print, "Making visible catalog...."
    obs_visibleone = get_bcc_obs_structure()
;    obs_visibleone=create_struct('id',0L, 'index', 0L, 'mag_u', 0.,$
;                                 'mag_g',0.,'mag_r',0.,'mag_i',0.,'mag_z',0.,$
;                                 'mag_y',0., 'magerr_g', 0., 'magerr_r', 0., $
;                                 'magerr_i', 0., 'magerr_z', 0., 'magerr_y', 0., $
;  			       'flux_g',0.,'flux_r',0.,'flux_i',0.,'flux_z',0.,'flux_y',0.,$
;  			       'ivar_g', 0., 'ivar_r', 0., 'ivar_i', 0., $
;  			       'ivar_z',0., 'ivar_y', 0., $
;                                 'ra',0.0,'dec',0.0,'epsilon', fltarr(2), 'size', 0., $
;                                 'zCarlos',0.0, 'ArborZ', 0., $
;                                 'ArborZ_err', 0., 'ANNZ', 0., 'ANNZ_err', 0.,$
;                                 'photoz_Gaussian', 0.)
    obs_visible = replicate(obs_visibleone, ng)

    print, "Filling the visible catalog..."
    obs_visible.id = obs_full[ii].id
    obs_visible.index = obs_full[ii].index
    obs_visible.ra = obs_full[ii].ra
    obs_visible.dec = obs_full[ii].dec
    obs_visible.mag_u = obs_full[ii].mag_u
    obs_visible.mag_g = obs_full[ii].omag[0]
    obs_visible.mag_r = obs_full[ii].omag[1]
    obs_visible.mag_i = obs_full[ii].omag[2]
    obs_visible.mag_z = obs_full[ii].omag[3]
    obs_visible.mag_y = obs_full[ii].omag[4]
    obs_visible.magerr_g = obs_full[ii].omagerr[0]
    obs_visible.magerr_r = obs_full[ii].omagerr[1]
    obs_visible.magerr_i = obs_full[ii].omagerr[2]
    obs_visible.magerr_z = obs_full[ii].omagerr[3]
    obs_visible.magerr_y = obs_full[ii].omagerr[4]
    obs_visible.flux_g = obs_full[ii].flux[0]
    obs_visible.flux_r = obs_full[ii].flux[1]
    obs_visible.flux_i = obs_full[ii].flux[2]
    obs_visible.flux_z = obs_full[ii].flux[3]
    obs_visible.flux_y = obs_full[ii].flux[4]
    obs_visible.ivar_g = obs_full[ii].ivar[0]
    obs_visible.ivar_g = obs_full[ii].ivar[0]
    obs_visible.ivar_g = obs_full[ii].ivar[0]
    obs_visible.ivar_g = obs_full[ii].ivar[0]
    obs_visible.epsilon = obs_full[ii].epsilon
    obs_visible.size = obs_full[ii].size
;   for i = 0L, 4 do begin
;     good = where(obs_full.omag[i] le mag_lim[i], comp=bad)
;     if (i eq 0 and bad[0] ge 0) then obs_visible[bad].mag_g = 99
;     if (i eq 1 and bad[0] ge 0) then obs_visible[bad].mag_r = 99
;     if (i eq 2 and bad[0] ge 0) then obs_visible[bad].mag_i = 99
;     if (i eq 3 and bad[0] ge 0) then obs_visible[bad].mag_z = 99
;     if (i eq 4 and bad[0] ge 0) then obs_visible[bad].mag_y = 99
;   endfor
    zerr=0.03*randomn(seed,N_ELEMENTS(obs_full))
    obs_full.photoz_Gaussian=obs_full.z+zerr*(1.0+obs_full.z)
    obs_visible.photoz_Gaussian = obs_full[ii].photoz_Gaussian
  
    ;;;define header information
    sxaddpar, hdr_t, 'EXTNAME', 'DES'
    sxaddpar, hdr_t, 'IND_G',0
    sxaddpar, hdr_t, 'IND_R',1
    sxaddpar, hdr_t, 'IND_I',2
    sxaddpar, hdr_t, 'IND_Z',3
    sxaddpar, hdr_t, 'IND_Y',4
    sxaddpar, hdr_t, 'LIMMAG_G', mag_lim5[0]
    sxaddpar, hdr_t, 'LIMMAG_R', mag_lim5[1]
    sxaddpar, hdr_t, 'LIMMAG_I', mag_lim5[2]
    sxaddpar, hdr_t, 'LIMMAG_Z', mag_lim5[3]
    sxaddpar, hdr_t, 'LIMMAG_Y', mag_lim5[4]
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
    sxaddpar, hdr_v, 'LIMMAG_G', mag_lim5[0]
    sxaddpar, hdr_v, 'LIMMAG_R', mag_lim5[1]
    sxaddpar, hdr_v, 'LIMMAG_I', mag_lim5[2]
    sxaddpar, hdr_v, 'LIMMAG_Z', mag_lim5[3]
    sxaddpar, hdr_v, 'LIMMAG_Y', mag_lim5[4]
    sxaddpar, hdr_v, 'AREA', 220.552
    sxaddpar, hdr_v, 'RA_MIN', 10.0
    sxaddpar, hdr_v, 'RA_MAX', 30.0
    sxaddpar, hdr_v, 'DEC_MIN', 35.0
    sxaddpar, hdr_v, 'DEC_MAX', 50.0

    sxaddpar, hdr_sdss, 'EXTNAME', 'SDSS'
    sxaddpar, hdr_sdss, 'IND_U',0
    sxaddpar, hdr_sdss, 'IND_G',1
    sxaddpar, hdr_sdss, 'IND_R',2
    sxaddpar, hdr_sdss, 'IND_I',3
    sxaddpar, hdr_sdss, 'IND_Z',4

    sxaddpar, hdr_vista, 'EXTNAME', 'Vista'
    sxaddpar, hdr_vista, 'IND_J',0
    sxaddpar, hdr_vista, 'IND_H',1
    sxaddpar, hdr_vista, 'IND_K',2

    sxaddpar, hdr_johnson, 'EXTNAME', 'Johnson'
    sxaddpar, hdr_johnson, 'IND_U',0
    sxaddpar, hdr_johnson, 'IND_B',1
    sxaddpar, hdr_johnson, 'IND_V',2
    sxaddpar, hdr_johnson, 'IND_R',3
    sxaddpar, hdr_johnson, 'IND_I',4

    sxaddpar, hdr_vista, 'EXTNAME', 'DEEP2'
    sxaddpar, hdr_vista, 'IND_B',0
    sxaddpar, hdr_vista, 'IND_R',1
    sxaddpar, hdr_vista, 'IND_I',2

    sxaddpar, hdr_flamex, 'EXTNAME', 'FLAMEX'
    sxaddpar, hdr_flamex, 'IND_Bw',0
    sxaddpar, hdr_flamex, 'IND_R',1
    sxaddpar, hdr_flamex, 'IND_I',2
    sxaddpar, hdr_flamex, 'IND_J',3
    sxaddpar, hdr_flamex, 'IND_H',4
    sxaddpar, hdr_flamex, 'IND_Ks',5
    sxaddpar, hdr_flamex, 'IND_3.6', 6
    sxaddpar, hdr_flamex, 'IND_4.5', 7

    sxaddpar, hdr_cfhtls, 'EXTNAME', 'CFHTLS'
    sxaddpar, hdr_cfhtls, 'IND_U', 0
    sxaddpar, hdr_cfhtls, 'IND_G', 0
    sxaddpar, hdr_cfhtls, 'IND_R', 0
    sxaddpar, hdr_cfhtls, 'IND_I', 0
    sxaddpar, hdr_cfhtls, 'IND_Z', 0

    sxaddpar, hdr_euclid, 'EXTNAME', 'EUCLID'
    sxaddpar, hdr_euclid, 'IND_VIS', 0
    sxaddpar, hdr_euclid, 'IND_Y', 1
    sxaddpar, hdr_euclid, 'IND_J', 2
    sxaddpar, hdr_euclid, 'IND_H', 3

    sxaddpar, hdr_irac, 'EXTNAME', 'IRAC'
    sxaddpar, hdr_irac, 'CHANNEL_1', 0
    sxaddpar, hdr_irac, 'CHANNEL_2', 1
    sxaddpar, hdr_irac, 'CHANNEL_3', 2
    sxaddpar, hdr_irac, 'CHANNEL_4', 3

    sxaddpar, hdr_wise, 'EXTNAME', 'WISE'
    sxaddpar, hdr_wise, 'W1', 0
    sxaddpar, hdr_wise, 'W2', 1
    sxaddpar, hdr_wise, 'W3', 2
    sxaddpar, hdr_wise, 'W4', 3

    sxaddpar, hdr_hsc, 'EXTNAME', 'HSC'
    sxaddpar, hdr_hsc, 'IND_G', 0
    sxaddpar, hdr_hsc, 'IND_R', 1
    sxaddpar, hdr_hsc, 'IND_I', 2
    sxaddpar, hdr_hsc, 'IND_Z', 3
    sxaddpar, hdr_hsc, 'IND_Y', 4

    sxaddpar, hdr_lsst, 'EXTNAME', 'LSST'
    sxaddpar, hdr_lsst, 'IND_U',0
    sxaddpar, hdr_lsst, 'IND_G',1
    sxaddpar, hdr_lsst, 'IND_R',2
    sxaddpar, hdr_lsst, 'IND_I',3
    sxaddpar, hdr_lsst, 'IND_Z',4
    sxaddpar, hdr_lsst, 'IND_Y2',5
    sxaddpar, hdr_lsst, 'IND_Y3',6
    sxaddpar, hdr_lsst, 'IND_Y4',7

    sxaddpar, hdr_wfirst, 'EXTNAME', 'WFIRST'
    sxaddpar, hdr_wfirst, 'IND_Y',0
    sxaddpar, hdr_wfirst, 'IND_J',1
    sxaddpar, hdr_wfirst, 'IND_H',2
    sxaddpar, hdr_wfirst, 'IND_K',3

    print, "Writing out first set of catalogs...."
    mwrfits, obs_full[ii], outpath+fbase+'_truth_no_photoz.'+healpix_num+'.fit', hdr_t, /create
    mwrfits, obs_visible, outpath+fbase+'_no_photoz.'+healpix_num+'.fit', hdr_v, /create
    mwrfits, sdss_mag[ii], outpath+fbase+'_sdss_mag.'+healpix_num+'.fit', hdr_sdss, /create
    mwrfits, vista_mag[ii], outpath+fbase+'_vista_mag.'+healpix_num+'.fit', hdr_vista, /create
    mwrfits, johnson_mag[ii], outpath+fbase+'_johnson_mag.'+healpix_num+'.fit', hdr_johnson, /create
    mwrfits, deep_mag[ii], outpath+fbase+'_deep_mag.'+healpix_num+'.fit', hdr_deep, /create
    mwrfits, flamex_mag[ii], outpath+fbase+'_flamex_mag.'+healpix_num+'.fit', hdr_flamex, /create
    mwrfits, cfhtls_mag[ii], outpath+fbase+'_cfhtls_mag.'+healpix_num+'.fit', hdr_cfhtls, /create
    mwrfits, euclid_mag[ii], outpath+fbase+'_euclid_mag.'+healpix_num+'.fit', hdr_euclid, /create
    mwrfits, irac_mag[ii], outpath+fbase+'_irac_mag.'+healpix_num+'.fit', hdr_irac, /create
    mwrfits, wise_mag[ii], outpath+fbase+'_wise_mag.'+healpix_num+'.fit', hdr_wise, /create
    mwrfits, hsc_mag[ii], outpath+fbase+'_hsc_mag.'+healpix_num+'.fit', hdr_hsc, /create
    mwrfits, lsst_mag[ii], outpath+fbase+'_lsst_mag.'+healpix_num+'.fit', hdr_lsst, /create
    mwrfits, wfirst_mag[ii], outpath+fbase+'_wfirst_mag.'+healpix_num+'.fit', hdr_wforst, /create
;  mwrfits, halos, outpath+fbase+'_halos.'+healpix_num+'.fit', /create

;  print, "Doing Photo-z's..."
;  add_bcc_zcarlos, path
  endif ;;;making sure there were legitimate galaxies to output


end



