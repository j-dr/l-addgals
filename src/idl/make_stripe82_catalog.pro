pro make_stripe82_catalog,path, hfile, mag_lim=mag_lim, no_lensing=no_lensing

  path += '/'

  g = mrdfits(path+'/galaxies_full.fit', 1)
  des = mrdfits(path+'/des_full.fit',1)
  vista = mrdfits(path+'/vista_full.fit',1)
  deep = mrdfits(path+'/deep_full.fit',1)
  johnson = mrdfits(path+'/johnson_full.fit',1)
  flamex = mrdfits(path+'/flamex_full.fit',1)
  sdss25 = mrdfits(path+'/sdss_full.fit',1)
  if NOT KEYWORD_SET(no_lensing) then begin
;    lg = mrdfits(path+'/from_joerg/galaxies_full_2.fit',1)
    lg = mrdfits(path+'/from_joerg/galaxies_full.fit',1)
    ldes = mrdfits(path+'/from_joerg/des_full.fit',1)
    lvista = mrdfits(path+'/from_joerg/vista_full.fit',1)
    ldeep = mrdfits(path+'/from_joerg/deep_full.fit',1)
    ljohnson = mrdfits(path+'/from_joerg/johnson_full.fit',1)
    lflamex = mrdfits(path+'/from_joerg/flamex_full.fit',1)
  endif

  if not KEYWORD_SET(mag_lim) then mag_lim = [22.0, 22.2, 22.2, 21.3, 20.5]+2.0

  ng=n_elements(g)
  if KEYWORD_SET(no_lensing) then ngnew = n_elements(g) else $
    ngnew = n_elements(lg)

  cooper_cat = mrdfits('/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/cooper/combined_dr6_cooper.fit',1)

  print, "Read the data, now making initial structures..."
  obs_fullone=create_struct('id',0L,'index', 0L, 'ecatid', 0L, 'coeffs', fltarr(5), $
                            'tmag',fltarr(5),'omag',fltarr(5), $
			    'flux', fltarr(5), 'ivar', fltarr(5), $
                            'omagerr',fltarr(5),'amag', fltarr(5), $
                            'ra',0.0,'dec',0.0,'z',0.0,$
                            'haloid',0L,'rhalo',0.0, 'm200', 0., 'ngals', 0L, $
                            'r200',0.0,'central',0, 'tra', 0., 'tdec', 0., $
                            'epsilon', fltarr(2), 'gamma1', 0., 'gamma2', 0., $
                            'kappa', 0., 'mu', 0., 'lmag', fltarr(5),$
                            'te', fltarr(2), 'tsize', 0., 'size', 0., $
                            'zCarlos',0.0, 'ArborZ', 0., $
                            'ArborZ_err', 0., 'ANNZ', 0., 'ANNZ_err', 0.,$
                            'photoz_Gaussian', 0., $
                            'px', 0., 'py', 0., 'pz', 0., $
                            'vx', 0., 'vy', 0., 'vz', 0.)
  des_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5), 'lmag', fltarr(5))
  vista_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5), 'lmag', fltarr(5))
  johnson_mag1 = create_struct('amag', fltarr(5), 'tmag', fltarr(5), 'lmag', fltarr(5))
  deep_mag1 = create_struct('amag', fltarr(3), 'tmag', fltarr(3), 'lmag', fltarr(3))
  flamex_mag1 = create_struct('amag', fltarr(8), 'tmag', fltarr(8), 'lmag', fltarr(8))
  sdss25_mag1 = create_struct('amag', fltarr(5))

  obs_full=replicate(obs_fullone,ngnew)
  des_mag = replicate(des_mag1, ngnew)
  vista_mag = replicate(vista_mag1, ngnew)
  johnson_mag = replicate(johnson_mag1, ngnew)
  deep_mag = replicate(deep_mag1, ngnew)
  flamex_mag = replicate(flamex_mag1, ngnew)
  sdss25_mag = replicate(sdss25_mag1, ngnew)

  ind = lindgen(ng)
  obs_full[ind].id=g.id
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
  obs_full[ind].tmag(0) = g.omag(0)
  obs_full[ind].tmag(1) = g.omag(1)
  obs_full[ind].tmag(2) = g.omag(2)
  obs_full[ind].tmag(3) = g.omag(3)
  obs_full[ind].tmag(4) = g.omag(4)
  obs_full[ind].amag(0) = g.amag(0)
  obs_full[ind].amag(1) = g.amag(1)
  obs_full[ind].amag(2) = g.amag(2)
  obs_full[ind].amag(3) = g.amag(3)
  obs_full[ind].amag(4) = g.amag(4)
  des_mag[ind].amag = des.amag
  des_mag[ind].tmag = des.omag
  vista_mag[ind].amag = vista.amag
  vista_mag[ind].tmag = vista.omag
  johnson_mag[ind].amag = johnson.amag
  johnson_mag[ind].tmag = johnson.omag
  deep_mag[ind].amag = deep.amag
  deep_mag[ind].tmag = deep.omag
  flamex_mag[ind].amag = flamex.amag
  flamex_mag[ind].tmag = flamex.omag
  sdss25_mag[ind].amag = sdss25.amag

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
       obs_full[i].tmag(0) = g[match].omag(0)
       obs_full[i].tmag(1) = g[match].omag(1)
       obs_full[i].tmag(2) = g[match].omag(2)
       obs_full[i].tmag(3) = g[match].omag(3)
       obs_full[i].tmag(4) = g[match].omag(4)
       obs_full[i].amag(0) = g[match].amag(0)
       obs_full[i].amag(1) = g[match].amag(1)
       obs_full[i].amag(2) = g[match].amag(2)
       obs_full[i].amag(3) = g[match].amag(3)
       obs_full[i].amag(4) = g[match].amag(4)
       vista_mag[i].amag = vista[match].amag
       vista_mag[i].tmag = vista[match].omag
       johnson_mag[i].amag = johnson[match].amag
       johnson_mag[i].tmag = johnson[match].omag
       deep_mag[i].amag = deep[match].amag
       deep_mag[i].tmag = deep[match].omag
       flamex_mag[i].amag = flamex[match].amag
       flamex_mag[i].tmag = flamex[match].omag
       des_mag[i].amag = des[match].amag
       des_mag[i].tmag = des[match].omag
       sdss25_mag[i].amag = sdss25[match].amag
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
    obs_full.lmag(0) = lg.lmag(0)
    obs_full.lmag(1) = lg.lmag(1)
    obs_full.lmag(2) = lg.lmag(2)
    obs_full.lmag(3) = lg.lmag(3)
    obs_full.lmag(4) = lg.lmag(4)
    sdss_mag.lmag = lg.lmag
    vista_mag.lmag = lvista.lmag
    johnson_mag.lmag = ljohnson.lmag
    deep_mag.lmag = ldeep.lmag
    flamex_mag.lmag = lflamex.lmag
  endif
  print, "Finished filling initial catalogs."
  
  ;;;delete some variables to save memory
  ldes = 0
  lg = 0
  lvista = 0
  ljohnson = 0
  ldeep = 0
  cooper_cat = 0

;  if KEYWORD_SET(no_lensing) then add_des_photometric_errors,obs_full, ka, obs_full2, maglim=mag_lim, /no_new_struct
;  if NOT KEYWORD_SET(no_lensing) then add_des_photometric_errors,obs_full, ka, obs_full2, maglim=mag_lim, /no_new_struct, /lmag
;  obs_full = obs_full2

;if(0) then begin
  max_galaxies_at_once = 10000
  galaxies_calculated = 0L
  ng = N_ELEMENTS(obs_full)
  while(galaxies_calculated lt ng) do begin
    n_to_do = min([max_galaxies_at_once, ng-galaxies_calculated])
    ind = galaxies_calculated + lindgen(n_to_do)
    if KEYWORD_SET(no_lensing) then add_sdss_photometric_errors,obs_full[ind], $
	obs_full2, maglim=mag_lim, /no_new_struct
    if NOT KEYWORD_SET(no_lensing) then add_sdss_photometric_errors,obs_full[ind], $
	obs_full2, maglim=mag_lim, /no_new_struct, /lmag
    obs_full[ind] = obs_full2[*]
    galaxies_calculated += n_to_do
    print, "Added errors for ", galaxies_calculated, " of ", ng, " galaxies."
  endwhile
;endif

  print, "Finished adding photometric errors."

  print, "Cutting the catalog..."
  ii = where(obs_full.omag[0] le mag_lim[0] or $
             obs_full.omag[1] le mag_lim[1] or $
             obs_full.omag[2] le mag_lim[2] or $
             obs_full.omag[3] le mag_lim[3] or $
             obs_full.omag[4] le mag_lim[4], count)
  obs_full = obs_full[ii]
  des_mag = des_mag[ii]
  vista_mag = vista_mag[ii]
  johnson_mag = johnson_mag[ii]
  deep_mag = deep_mag[ii]
  flamex_mag = flamex_mag[ii]
  sdss25_mag = sdss25_mag[ii]
  ng = count
;  for i = 0L, 4 do begin
;     good = where(obs_full.omag[i] le mag_lim[i], comp=bad)
;     if (bad[0] ge 0) then obs_full[bad].omag[i] = 99
;  endfor
  
  print, "Adding HOD info...."
  fix_hod, halos, obs_full, hfile=hfile, /ascii, boxsize=3000., llbins=128, /non_sham, /addgals

  zerr=0.05*randomn(seed,ng)
  obs_full.photoz_Gaussian=obs_full.z+zerr*(1.0+obs_full.z)
;;;  obs_visible.photoz_Gaussian = obs_full.photoz_Gaussian

  ;;;define header information
  sxaddpar, hdr_t, 'EXTNAME', 'SDSS'
  sxaddpar, hdr_t, 'IND_U',0
  sxaddpar, hdr_t, 'IND_G',1
  sxaddpar, hdr_t, 'IND_R',2
  sxaddpar, hdr_t, 'IND_I',3
  sxaddpar, hdr_t, 'IND_Z',4
  sxaddpar, hdr_t, 'LIMMAG_U', mag_lim[0]
  sxaddpar, hdr_t, 'LIMMAG_G', mag_lim[1]
  sxaddpar, hdr_t, 'LIMMAG_R', mag_lim[2]
  sxaddpar, hdr_t, 'LIMMAG_I', mag_lim[3]
  sxaddpar, hdr_t, 'LIMMAG_Z', mag_lim[4]
  sxaddpar, hdr_t, 'AREA', 220.552
  sxaddpar, hdr_t, 'RA_MIN', 10.0
  sxaddpar, hdr_t, 'RA_MAX', 30.0
  sxaddpar, hdr_t, 'DEC_MIN', 35.0
  sxaddpar, hdr_t, 'DEC_MAX', 50.0
  sxaddpar, hdr_t, 'Z_LO', 0.0
  sxaddpar, hdr_t, 'Z_HI', 1.33


  sxaddpar, hdr_des, 'EXTNAME', 'DES'
  sxaddpar, hdr_des, 'IND_G',0
  sxaddpar, hdr_des, 'IND_R',1
  sxaddpar, hdr_des, 'IND_I',2
  sxaddpar, hdr_des, 'IND_Z',3
  sxaddpar, hdr_des, 'IND_Y',4

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

  print, "Writing out first set of catalogs...."
  mwrfits, obs_full, path+'Stripe82_Mock_no_photoz.fit', hdr_t, /create
  mwrfits, dees_mag, path+'Stripe82_Mock_des_mag.fit', hdr_des, /create
  mwrfits, vista_mag, path+'Stripe82_Mock_vista_mag.fit', hdr_vista, /create
  mwrfits, johnson_mag, path+'Stripe82_Mock_johnson_mag.fit', hdr_johnson, /create
  mwrfits, deep_mag, path+'Stripe82_Mock_deep_mag.fit', hdr_deep, /create
  mwrfits, flamex_mag, path+'Stripe82_Mock_flamex_mag.fit', hdr_deep, /create
  mwrfits, sdss25_mag, path+'Stripe82_Mock_sdss25_mag.fit', hdr_deep, /create
  mwrfits, halos, path+'Stripe82_Mock_halos.fit', /create

asdf
  print, "Getting photo-z training set..."
;  ii = get_training_set(obs_full, random=0, nrandom = 3000, /stripe82, deep=deep, /vvds, /zcosmos)
  ii = get_training_set(obs_full, random=0, nrandom = 3000, /stripe82, deep=deep_mag, /vvds, /zcosmos, sdss=sdss_mag)
  srt = sort(obs_full[ii].ra)
  validation = ii[srt[0:9999]]
  training = ii[srt[10000:N_ELEMENTS(ii)-1]]
  mwrfits, obs_full[validation], path+'/DES_Mock_validation_set.fit', /create
  mwrfits, obs_full[training], path+'/DES_Mock_training_set.fit', /create

end

  print, "Doing zCarlos..."
  make_carlos_photoz, path+'/DES_Mock_training_set.fit', $
                      path+'/DES_Mock_Baseline_truth.fit', $
                      ofile = path+'/DES_Mock_zCarlos.fit'
  
  zCarlos = mrdfits(path+'/DES_Mock_zCarlos.fit',1)
  obs_full.zCarlos = zCarlos.zCarlos_photoz
  obs_visible.zCarlos = zCarlos.zCarlos_photoz

  print, "Writing sub-regions..."
  ;;make sub regions that are just 2 degrees wide
  ramin = 10
  dra = 2
  for i = 0, 9 do begin
     tramin = ramin + dra*i
     tramax = tramin + dra
     ii = where(obs_full.ra ge tramin and obs_full.ra lt tramax)

     sxdelpar, hdr_t, 'RA_MIN'
     sxdelpar, hdr_t, 'RA_MAX'
;     delvar, hdr_t
;     sxaddpar, hdr_t, 'EXTNAME', 'DES'
;     sxaddpar, hdr_t, 'IND_G',0
;     sxaddpar, hdr_t, 'IND_R',1
;     sxaddpar, hdr_t, 'IND_I',2
;     sxaddpar, hdr_t, 'IND_Z',3
;     sxaddpar, hdr_t, 'IND_Y',4
;     sxaddpar, hdr_t, 'LIMMAG_G', 24.6
;     sxaddpar, hdr_t, 'LIMMAG_R', 24.1
;     sxaddpar, hdr_t, 'LIMMAG_I', 24.4
;     sxaddpar, hdr_t, 'LIMMAG_Z', 23.8
;     sxaddpar, hdr_t, 'LIMMAG_Y', 21.3
;     sxaddpar, hdr_t, 'AREA', 22.0552
     sxaddpar, hdr_t, 'RA_MIN', tramin
     sxaddpar, hdr_t, 'RA_MAX', tramax
;     sxaddpar, hdr_t, 'DEC_MIN', 35.0
;     sxaddpar, hdr_t, 'DEC_MAX', 50.0
;     sxaddpar, hdr_t, 'Z_LO', 0.0
;     sxaddpar, hdr_t, 'Z_HI', 1.33
     
     sxdelpar, hdr_v, 'RA_MIN'
     sxdelpar, hdr_v, 'RA_MAX'
;     delvar, hdr_v
;     sxaddpar, hdr_v, 'EXTNAME', 'DES'
;     sxaddpar, hdr_v, 'IND_G',0
;     sxaddpar, hdr_v, 'IND_R',1
;     sxaddpar, hdr_v, 'IND_I',2
;     sxaddpar, hdr_v, 'IND_Z',3
;     sxaddpar, hdr_v, 'IND_Y',4
;     sxaddpar, hdr_v, 'LIMMAG_G', 24.6
;     sxaddpar, hdr_v, 'LIMMAG_R', 24.1
;     sxaddpar, hdr_v, 'LIMMAG_I', 24.4
;     sxaddpar, hdr_v, 'LIMMAG_Z', 23.8
;     sxaddpar, hdr_v, 'LIMMAG_Y', 21.3
;     sxaddpar, hdr_v, 'AREA', 22.0552
     sxaddpar, hdr_v, 'RA_MIN', tramin
     sxaddpar, hdr_v, 'RA_MAX', tramax
;     sxaddpar, hdr_v, 'DEC_MIN', 35.0
;     sxaddpar, hdr_v, 'DEC_MAX', 50.0

     mwrfits, obs_full[ii], path+'/DES_Mock_Baseline_truth.'+strcompress(string(i),/remove_all)+'.fit', hdr_t, /create
     mwrfits, obs_visible[ii], path+'/DES_Mock_Baseline.'+strcompress(string(i),/remove_all)+'.fit', hdr_v, /create
     mwrfits, zCarlos[ii], path+'/DES_Mock_zCarlos.'+strcompress(string(i),/remove_all)+'.fit', /create
  endfor

  mwrfits, obs_full, path+'DES_Mock_Baseline_truth.fit', truth_hdr, /create
  mwrfits, obs_visible, path+'DES_Mock_Baseline.fit', obs_hdr, /create



end



