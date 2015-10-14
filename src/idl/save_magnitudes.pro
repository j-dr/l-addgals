pro save_magnitudes, name, gals, dim=dim, $
                    des=des, $
                    vista=vista, deep=deep,johnson=johnson, ubv_deep=ubv_deep,$
                    alhambra=alhambra,$
                    flamex=flamex, cfhtls=cfhtls, hsc=hsc, lsst=lsst, $
		    euclid=euclid, wfirst=wfirst, twomass=twomass, irac=irac, wise=wise

;;;write out the new magnitude
IF (keyword_set(des)) then begin
  desname=name+'_des'
  IF keyword_set(dim) then desname += '_dim'
  desname += '.fit'
  gals_des1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_des = replicate(gals_des1,N_ELEMENTS(gals))
  gals_des.omag = gals.odesmag
  gals_des.amag = gals.adesmag
  mwrfits, gals_des, desname, /create
ENDIF
IF keyword_set(vista) then begin
  vistaname=name+'_vista'
  IF keyword_set(dim) then vistaname += '_dim'
  vistaname += '.fit'
  gals_vista1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_vista = replicate(gals_vista1,N_ELEMENTS(gals))
  gals_vista.omag = gals.ovistamag
  gals_vista.amag = gals.avistamag
  mwrfits, gals_vista, vistaname, /create
ENDIF
IF keyword_set(deep) then begin
  deepname=name+'_deep'
  IF keyword_set(dim) then deepname += '_dim'
  deepname += '.fit'
  gals_deep1 = create_struct('omag',fltarr(3), 'amag',fltarr(3))
  gals_deep = replicate(gals_deep1,N_ELEMENTS(gals))
  gals_deep.omag = gals.odeepmag
  gals_deep.amag = gals.adeepmag
  mwrfits, gals_deep, deepname, /create
ENDIF
IF keyword_set(johnson) then begin
  johnsonname=name+'_johnson'
  IF keyword_set(dim) then johnsonname += '_dim'
  johnsonname += '.fit'
  gals_johnson1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_johnson = replicate(gals_johnson1,N_ELEMENTS(gals))
  gals_johnson.omag = gals.ojohnsonmag
  gals_johnson.amag = gals.ajohnsonmag
  mwrfits, gals_johnson, johnsonname, /create
ENDIF
IF keyword_set(ubv_deep) then begin
  ubv_deepname=name+'_ubv_deep'
  IF keyword_set(dim) then ubv_deepname += '_dim'
  ubv_deepname += '.fit'
  gals_ubv_deep1 = create_struct('omag',fltarr(3), 'amag',fltarr(3))
  gals_ubv_deep = replicate(gals_ubv_deep1,N_ELEMENTS(gals))
  gals_ubv_deep.omag = gals.oubv_deepmag
  gals_ubv_deep.amag = gals.aubv_deepmag
  mwrfits, gals_ubv_deep, ubv_deepname, /create
ENDIF
IF keyword_set(alhambra) then begin
  alhambraname=name+'_alhambra'
  IF keyword_set(dim) then alhambraname += '_dim'
  alhambraname += '.fit'
  gals_alhambra1 = create_struct('omag',fltarr(20), 'amag',fltarr(20))
  gals_alhambra = replicate(gals_alhambra1,N_ELEMENTS(gals))
  gals_alhambra.omag = gals.oalhambramag
  gals_alhambra.amag = gals.aalhambramag
  mwrfits, gals_alhambra, alhambraname, /create
ENDIF
IF keyword_set(flamex) then begin
  flamexname=name+'_flamex'
  IF keyword_set(dim) then flamexname += '_dim'
  flamexname += '.fit'
  gals_flamex1 = create_struct('omag',fltarr(8), 'amag',fltarr(8))
  gals_flamex = replicate(gals_flamex1,N_ELEMENTS(gals))
  gals_flamex.omag = gals.oflamexmag
  gals_flamex.amag = gals.aflamexmag
  mwrfits, gals_flamex, flamexname, /create
ENDIF
IF keyword_set(cfhtls) then begin
  cfhtlsname=name+'_cfhtls'
  IF keyword_set(dim) then cfhtlsname += '_dim'
  cfhtlsname += '.fit'
  gals_cfhtls1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_cfhtls = replicate(gals_cfhtls1,N_ELEMENTS(gals))
  gals_cfhtls.omag = gals.ocfhtlsmag
  gals_cfhtls.amag = gals.acfhtlsmag
  mwrfits, gals_cfhtls, cfhtlsname, /create
ENDIF
IF keyword_set(hsc) then begin
  hscname=name+'_hsc'
  IF keyword_set(dim) then hscname += '_dim'
  hscname += '.fit'
  gals_hsc1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_hsc = replicate(gals_hsc1,N_ELEMENTS(gals))
  gals_hsc.omag = gals.ohscmag
  gals_hsc.amag = gals.ahscmag
  mwrfits, gals_hsc, hscname, /create
ENDIF
IF keyword_set(lsst) then begin
  lsstname=name+'_lsst'
  IF keyword_set(dim) then lsstname += '_dim'
  lsstname += '.fit'
  gals_lsst1 = create_struct('omag',fltarr(8), 'amag',fltarr(8))
  gals_lsst = replicate(gals_lsst1,N_ELEMENTS(gals))
  gals_lsst.omag = gals.olsstmag
  gals_lsst.amag = gals.alsstmag
  mwrfits, gals_lsst, lsstname, /create
ENDIF
IF keyword_set(euclid) then begin
  euclidname=name+'_euclid'
  IF keyword_set(dim) then euclidname += '_dim'
  euclidname += '.fit'
  gals_euclid1 = create_struct('omag',fltarr(4), 'amag',fltarr(4))
  gals_euclid = replicate(gals_euclid1,N_ELEMENTS(gals))
  gals_euclid.omag = gals.oeuclidmag
  gals_euclid.amag = gals.aeuclidmag
  mwrfits, gals_euclid, euclidname, /create
ENDIF
IF keyword_set(irac) then begin
  iracname=name+'_irac'
  IF keyword_set(dim) then iracname += '_dim'
  iracname += '.fit'
  gals_irac1 = create_struct('omag',fltarr(4), 'amag',fltarr(4))
  gals_irac = replicate(gals_irac1,N_ELEMENTS(gals))
  gals_irac.omag = gals.oiracmag
  gals_irac.amag = gals.airacmag
  mwrfits, gals_irac, iracname, /create
ENDIF
IF keyword_set(wise) then begin
  wisename=name+'_wise'
  IF keyword_set(dim) then wisename += '_dim'
  wisename += '.fit'
  gals_wise1 = create_struct('omag',fltarr(4), 'amag',fltarr(4))
  gals_wise = replicate(gals_wise1,N_ELEMENTS(gals))
  gals_wise.omag = gals.owisemag
  gals_wise.amag = gals.awisemag
  mwrfits, gals_wise, wisename, /create
ENDIF
IF keyword_set(wfirst) then begin
  wfirstname=name+'_wfirst'
  IF keyword_set(dim) then wfirstname += '_dim'
  wfirstname += '.fit'
  gals_wfirst1 = create_struct('omag',fltarr(4), 'amag',fltarr(4))
  gals_wfirst = replicate(gals_wfirst1,N_ELEMENTS(gals))
  gals_wfirst.omag = gals.owfirstmag
  gals_wfirst.amag = gals.awfirstmag
  mwrfits, gals_wfirst, wfirstname, /create
ENDIF

print, "Finished saving the files."

end

