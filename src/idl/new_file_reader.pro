pro new_file_reader, str, pathname=pathname,des=des, vista=vista, deep=deep,$
                     johnson=johnson, ubv_deep=ubv_deep,dim=dim, $
                     alhambra=alhambra, flamex=flamex, cfhtls=cfhtls,$
		     hsc=hsc, lsst=lsst, euclid=euclid, wfirst=wfirst, $
                     twomass=twomass, irac=irac, wise=wise, candels=candels,$
                     galz=galz, cut=cut
;if not keyword_set(pathname) then pathname='./'

if (cut eq 0) then begin
   fname1=pathname+'gal_ginfo.dat'
endif else begin
   fname1=pathname+'gal_ginfo.dat.cut'
endelse
;   fname2=pathname+'gal_pinfo.dat'
;   fname3=pathname+'gal_hinfo.dat'

   spawn,'wc '+fname1,result
   rarray=str_sep(strtrim(result(0),2),' ')
   nobj=long(rarray(0))
   print, nobj
   s1=create_struct('id',0L,'ecatid', 0L, 'omag',fltarr(5),$
		    'amag',fltarr(5),'z',0.0,'ra',0.0,'dec',0.0,$
                    'px',0.0,'py',0.0,'pz',0.0,'vx',0.0,'vy',0.0,'vz',0.0,$
                    'edge', 0, 'siglos', 0.0, $
                    'halopx', 0.0, 'halopy', 0.0, 'halopz', 0.0, $
                    'halovx', 0.0, 'halovy', 0.0, 'halovz', 0.0, $
                    'haloid',0L,'m200',0.0,'ngals',0,'r200',0.0,$
		    'rhalo',0.0,'halora',0.0,'halodec',0.0,'haloz',0.0, $
   ;'d8', 0.0, 'dg', 0.0, 
                    'central',0, 'mstar', 0., 'd8', 0., 'd8_drawn', 0., 'nndist', 0., 'nnpercent', 0.)
print, 'struct created'
IF keyword_set(des) THEN $
  s1=create_struct(s1,'odesmag',fltarr(6),'adesmag',fltarr(6))
IF keyword_set(vista) THEN $
  s1=create_struct(s1,'ovistamag',fltarr(5),'avistamag',fltarr(5))
IF keyword_set(deep) THEN $
  s1=create_struct(s1,'odeepmag',fltarr(3),'adeepmag',fltarr(3))
IF keyword_set(johnson) THEN $
  s1=create_struct(s1,'ojohnsonmag',fltarr(5),'ajohnsonmag',fltarr(5))
IF keyword_set(ubv_deep) THEN $
  s1=create_struct(s1,'oubv_deepmag',fltarr(3),'aubv_deepmag',fltarr(3))
IF keyword_set(alhambra) THEN $
  s1=create_struct(s1,'oalhambramag',fltarr(20),'aalhambramag',fltarr(20))
IF keyword_set(flamex) THEN $
  s1=create_struct(s1,'oflamexmag',fltarr(8),'aflamexmag',fltarr(8))
IF keyword_set(cfhtls) THEN $
  s1=create_struct(s1,'ocfhtlsmag',fltarr(5),'acfhtlsmag',fltarr(5))
IF keyword_set(hsc) THEN $
  s1=create_struct(s1,'ohscmag',fltarr(5),'ahscmag',fltarr(5))
IF keyword_set(lsst) THEN $
  s1=create_struct(s1,'olsstmag',fltarr(8),'alsstmag',fltarr(8))
IF keyword_set(euclid) THEN $
  s1=create_struct(s1,'oeuclidmag',fltarr(4),'aeuclidmag',fltarr(4))
IF keyword_set(irac) THEN $
  s1=create_struct(s1,'oiracmag',fltarr(4),'airacmag',fltarr(4))
IF keyword_set(wise) THEN $
  s1=create_struct(s1,'owisemag',fltarr(4),'awisemag',fltarr(4))
IF keyword_set(wfirst) THEN $
  s1=create_struct(s1,'owfirstmag',fltarr(4),'awfirstmag',fltarr(4))
IF keyword_set(twomass) THEN $
  s1 = create_struct(s1, 'otwomassmag', fltarr(3), 'atwomassmag', fltarr(3))
IF keyword_set(candels) THEN $
  s1 = create_struct(s1, 'ocandelsmag', fltarr(10), 'acandelsmag', fltarr(10))
IF keyword_set(galz) THEN $
  s1 = create_struct(s1, 'gz', 0.0, 'zbin', 0L, 'gzcentral', 0L)

   str=replicate(s1,nobj)
   str.id = lindgen(nobj)

   file1 = pathname+'gal_hinfo'
   file2 = pathname+'gal_pinfo'
   file3 = pathname+'gal_zinfo'
if (dim eq 1) then begin
   file1 += '_dim'
   file2 += '_dim'
   file3 += '_dim'
endif
   file1 += '.dat'
   file2 += '.dat'
   file3 += '.dat'
if (cut eq 1) then begin
   file1 += '.cut'
   file2 += '.cut'
   file3 += '.cut'
endif

   rdfloat, file1, haloid, m200, ngals, r200, rhalo, siglos, $
            halopx, halopy, halopz, halovx, halovy, halovz, $
            halora, halodec, haloz
   print, 'read gal hinfo'
   rdfloat, file2, px, py, pz, vx, vy, vz
if keyword_set(galz) then begin
   print, 'read gal zinfo'
   rdfloat, file3, gz, zbin, gzcentral
endif

   print, 'read gal pinfo'

   rhalo_old = rhalo
   ng = N_ELEMENTS(rhalo)
   for i = 0L, ng - 1 do begin
     rhalo(i) = sqrt((px(i)-halopx(i))^2 + (py(i)-halopy(i))^2 + (pz(i)-halopz(i))^2)
   endfor

;   rdfloat, pathname+'gal_dinfo.dat', d
;   rdfloat, pathname+'dmdgmr.dat', d8, dg, mr
if (cut eq 0) then begin
   rdfloat, pathname+'gal_rnninfo.dat', nndist, nnpercent
   print, 'rnn info'
   rdfloat, pathname+'gal_dinfo.dat', td8, d8
   print, 'dinfo'

   rdfloat, pathname+'gal_ginfo.dat', u, g, r, i, z, mu, mg, mr, mi, mz, ra, dec, zz, central, ecatid
   print, 'ginfo'
endif else begin
   rdfloat, pathname+'gal_rnninfo.dat.cut', nndist, nnpercent
   print, 'rnn info'
   rdfloat, pathname+'gal_dinfo.dat.cut', td8, d8
   print, 'dinfo'

   rdfloat, pathname+'gal_ginfo.dat.cut', u, g, r, i, z, mu, mg, mr, mi, mz, ra, dec, zz, central, ecatid
   print, 'ginfo'
endelse
IF keyword_set(des) THEN begin
;  rdfloat, pathname+'gal_ginfo_des.dat', des_y, des_z, des_g, des_i, des_r, des_z2,$
;           mdes_y, mdes_z, mdes_g, mdes_i, mdes_r, mdes_z2
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_des.dat', des_u, des_g, des_r, des_i, des_z, des_y, $
               mdes_u, mdes_g, mdes_r, mdes_i, mdes_z, mdes_y
   endif else begin
      rdfloat, pathname+'gal_ginfo_des.dat.cut', des_u, des_g, des_r, des_i, des_z, des_y, $
               mdes_u, mdes_g, mdes_r, mdes_i, mdes_z, mdes_y
   endelse
   print, 'des'
endif
IF keyword_set(vista) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_vista.dat',wfcam_j,wfcam_h,wfcam_k,wfcam_y,wfcam_z,$
               mwfcam_j,mwfcam_h,mwfcam_k,mwfcam_y,mwfcam_z
   endif else begin
      rdfloat, pathname+'gal_ginfo_vista.dat.cut',wfcam_j,wfcam_h,wfcam_k,wfcam_y,wfcam_z,$
               mwfcam_j,mwfcam_h,mwfcam_k,mwfcam_y,mwfcam_z
   endelse
   print, 'vista'
endif
IF keyword_set(deep) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_deep.dat',deep_b,deep_r,deep_i,$
               mdeep_b, mdeep_r, mdeep_i
   endif else begin
      rdfloat, pathname+'gal_ginfo_deep.dat.cut',deep_b,deep_r,deep_i,$
               mdeep_b, mdeep_r, mdeep_i
   endelse
   print, 'deep'
endif
IF keyword_set(johnson) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_johnson.dat',johnson_U,johnson_B,johnson_V,$
               johnson_R,johnson_I,$
               mjohnson_U, mjohnson_B, mjohnson_V, mjohnson_R, mjohnson_I
   endif else begin
      rdfloat, pathname+'gal_ginfo_johnson.dat.cut',johnson_U,johnson_B,johnson_V,$
               johnson_R,johnson_I,$
               mjohnson_U, mjohnson_B, mjohnson_V, mjohnson_R, mjohnson_I
   endelse
   print, 'johnson'
endif
IF keyword_set(ubv_deep) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_ubv_deep.dat',ubv_deep_U,ubv_deep_B,ubv_deep_V,$
               mubv_deep_U, mubv_deep_B, mubv_deep_V
   endif else begin
      rdfloat, pathname+'gal_ginfo_ubv_deep.dat.cut',ubv_deep_U,ubv_deep_B,ubv_deep_V,$
               mubv_deep_U, mubv_deep_B, mubv_deep_V
   endelse   
   print, 'uvb_deep'
endif   
IF keyword_set(alhambra) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_alhambra.dat',filter1,filter2,filter3,filter4,filter5,filter6,$
               filter7,filter8,filter9,filter10,filter11,filter12,filter13,filter14,filter15,$
               filter16,filter17,filter18,filter19,filter20,$
               mfilter1,mfilter2,mfilter3,mfilter4,mfilter5,mfilter6,$
               mfilter7,mfilter8,mfilter9,mfilter10,mfilter11,mfilter12,mfilter13,mfilter14,mfilter15,$
               mfilter16,mfilter17,mfilter18,mfilter19,mfilter20
   endif else begin
      rdfloat, pathname+'gal_ginfo_alhambra.dat.cut',filter1,filter2,filter3,filter4,filter5,filter6,$
               filter7,filter8,filter9,filter10,filter11,filter12,filter13,filter14,filter15,$
               filter16,filter17,filter18,filter19,filter20,$
               mfilter1,mfilter2,mfilter3,mfilter4,mfilter5,mfilter6,$
               mfilter7,mfilter8,mfilter9,mfilter10,mfilter11,mfilter12,mfilter13,mfilter14,mfilter15,$
               mfilter16,mfilter17,mfilter18,mfilter19,mfilter20
   endelse
   print, 'alhambra'
endif
IF keyword_set(flamex) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_flamex.dat',flamex_Bw,flamex_R,flamex_I,$
               flamex_J, flamex_H, flamex_Ks, flamex_36, flamex_45, $
               mflamex_Bw, mflamex_R, mflamex_I, $
               mflamex_J, mflamex_H, mflamex_Ks, mflamex_36, mflamex_45
   endif else begin
      rdfloat, pathname+'gal_ginfo_flamex.dat.cut',flamex_Bw,flamex_R,flamex_I,$
               flamex_J, flamex_H, flamex_Ks, flamex_36, flamex_45, $
               mflamex_Bw, mflamex_R, mflamex_I, $
               mflamex_J, mflamex_H, mflamex_Ks, mflamex_36, mflamex_45
   endelse
   print, 'flamex'
endif
IF keyword_set(cfhtls) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_cfhtls.dat', cfhtls_u, cfhtls_g, cfhtls_r,$
               cfhtls_i, cfhtls_z, $
               mcfhtls_u, mcfhtls_g, mcfhtls_r, $
               mcfhtls_i, mcfhtls_z
   endif else begin
      rdfloat, pathname+'gal_ginfo_cfhtls.dat.cut', cfhtls_u, cfhtls_g, cfhtls_r,$
               cfhtls_i, cfhtls_z, $
               mcfhtls_u, mcfhtls_g, mcfhtls_r, $
               mcfhtls_i, mcfhtls_z
   endelse
   print, 'cfhtls'
endif
IF keyword_set(hsc) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_hsc.dat', hsc_g, hsc_r,$
               hsc_i, hsc_z, hsc_y, $
               mhsc_g, mhsc_r, $
               mhsc_i, mhsc_z, mhsc_y
   endif else begin
      rdfloat, pathname+'gal_ginfo_hsc.dat.cut', hsc_g, hsc_r,$
               hsc_i, hsc_z, hsc_y, $
               mhsc_g, mhsc_r, $
               mhsc_i, mhsc_z, mhsc_y
   endelse 
   print, 'hsc'
endif
IF keyword_set(lsst) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_lsst.dat', lsst_u, lsst_g, lsst_r,$
               lsst_i, lsst_z, lsst_y2, lsst_y3, lsst_y4, $
               mlsst_u, mlsst_g, mlsst_r, $
               mlsst_i, mlsst_z, mlsst_y2, mlsst_y3, mlsst_y4
   endif else begin
      rdfloat, pathname+'gal_ginfo_lsst.dat.cut', lsst_u, lsst_g, lsst_r,$
               lsst_i, lsst_z, lsst_y2, lsst_y3, lsst_y4, $
               mlsst_u, mlsst_g, mlsst_r, $
               mlsst_i, mlsst_z, mlsst_y2, mlsst_y3, mlsst_y4
   endelse
   print, 'lsst'
endif
IF keyword_set(euclid) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_euclid.dat', euclid_vis, euclid_y, euclid_j,$
               euclid_h, $
               meuclid_vis, meuclid_y, meuclid_j, $
               meuclid_h
   endif else begin
      rdfloat, pathname+'gal_ginfo_euclid.dat.cut', euclid_vis, euclid_y, euclid_j,$
               euclid_h, $
               meuclid_vis, meuclid_y, meuclid_j, $
               meuclid_h
   endelse
   print, 'euclid'
endif
IF keyword_set(irac) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_irac.dat', irac_1, irac_2, irac_3,$
               irac_4, $
               mirac_1, mirac_2, mirac_3, $
               mirac_4
   endif else begin
      rdfloat, pathname+'gal_ginfo_irac.dat.cut', irac_1, irac_2, irac_3,$
               irac_4, $
               mirac_1, mirac_2, mirac_3, $
               mirac_4
   endelse
   print, 'irac'
endif
IF keyword_set(wise) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_wise.dat', wise_1, wise_2, wise_3,$
               wise_4, $
               mwise_1, mwise_2, mwise_3, $
               mwise_4
   endif else begin
      rdfloat, pathname+'gal_ginfo_wise.dat.cut', wise_1, wise_2, wise_3,$
               wise_4, $
               mwise_1, mwise_2, mwise_3, $
               mwise_4
   endelse
   print, 'wise'
endif
IF keyword_set(wfirst) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_wfirst.dat', wfirst_y, wfirst_j, wfirst_h,$
               wfirst_k, $
               mwfirst_y, mwfirst_j, mwfirst_h, $
               mwfirst_k
   endif else begin
      rdfloat, pathname+'gal_ginfo_wfirst.dat.cut', wfirst_y, wfirst_j, wfirst_h,$
               wfirst_k, $
               mwfirst_y, mwfirst_j, mwfirst_h, $
               mwfirst_k
   endelse
   print, 'wfirst'
endif
IF keyword_set(twomass) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_twomass.dat', twomass_j, twomass_h,twomass_ks, $
               mtwomass_j, mtwomass_h, mtwomass_ks
   endif else begin
      rdfloat, pathname+'gal_ginfo_twomass.dat.cut', twomass_j, twomass_h,twomass_ks, $
               mtwomass_j, mtwomass_h, mtwomass_ks
   endelse
   print, 'twomass'
endif
IF keyword_set(candels) THEN begin
   if (cut eq 0) then begin
      rdfloat, pathname+'gal_ginfo_candels.dat', candels_0, candels_1, candels_2, $
               candels_3, candels_4, candels_5, $
               candels_6, candels_7, candels_8, $
               candels_9, candels_10, $
               mcandels_0, mcandels_1, mcandels_2, $
               mcandels_3, mcandels_4, mcandels_5, $
               mcandels_6, mcandels_7, mcandels_8, $
               mcandels_9, mcandels_10
   endif else begin
      rdfloat, pathname+'gal_ginfo_candels.dat.cut', candels_0, candels_1, candels_2, $
               candels_3, candels_4, candels_5, $
               candels_6, candels_7, candels_8, $
               candels_9,$
               mcandels_0, mcandels_1, mcandels_2, $
               mcandels_3, mcandels_4, mcandels_5, $
               mcandels_6, mcandels_7, mcandels_8, $
               mcandels_9
   endelse
   print, 'candels'
endif



IF keyword_set(declimit) THEN BEGIN
    ii = where(dec LT 5)
    haloid = haloid[ii]
    m200 = m200[ii]
    ngals = ngals[ii]
    r200 = r200[ii]
    rhalo = rhalo[ii]
    siglos = siglos[ii]
    halopx = halopx[ii]
    halopy = halopy[ii]
    halopz = halopz[ii]

    halovx = halovx[ii]
    halovy = halovy[ii]
    halovz = halovz[ii]

    px = px[ii]
    py = py[ii]
    pz = pz[ii]
    vx = vx[ii]
    vy = vy[ii]
    vz = vz[ii]

    id = id[ii]
    td8 = td8[ii]
    d8 = d8[ii]
    dg = dg[ii]
    u = u[ii]
    g = g[ii]
    r = r[ii]
    i = i[ii]
    z = z[ii]
    mu = mu[ii]
    mg = mg[ii]
    mr = mr[ii]
    mi = mi[ii]
    mz = mz[ii]
    ra = ra[ii]
    dec = dec[ii]
    zz = zz[ii]
    central = central[ii]
    ecatid = ecatid[ii]
    nndist = nndist[ii]
    nnpercent = nnpercent[ii]
    print, 'declimit'
IF keyword_set(des) then begin
    des_u = des_u[ii] 
    des_y = des_y[ii] 
    des_z = des_z[ii] 
    des_g = des_g[ii] 
    des_i = des_i[ii] 
    des_r = des_r[ii] 
;    des_z2 = des_z2[ii]
    mdes_u = mdes_u[ii] 
    mdes_y = mdes_y[ii] 
    mdes_z = mdes_z[ii] 
    mdes_g = mdes_g[ii] 
    mdes_i = mdes_i[ii] 
    mdes_r = mdes_r[ii] 
;    mdes_z2 = mdes_z2[ii]
    print, 'des'
ENDIF
IF keyword_set(vista) then begin
    wfcam_j = wfcam_j[ii]
    wfcam_h = wfcam_h[ii]
    wfcam_k = wfcam_k[ii]
    wfcam_y = wfcam_y[ii]
    wfcam_z = wfcam_z[ii]
    mwfcam_j = mwfcam_j[ii]
    mwfcam_h = mwfcam_h[ii]
    mwfcam_k = mwfcam_k[ii]
    mwfcam_y = mwfcam_y[ii]
    mwfcam_z = mwfcam_z[ii]
    print, 'vista'
ENDIF
IF keyword_set(deep) then begin
    deep_b = deep_b[ii]
    deep_r = deep_r[ii]
    deep_i = deep_i[ii]
    mdeep_b = mdeep_b[ii]
    mdeep_r = mdeep_r[ii]
    mdeep_i = mdeep_i[ii]
    print, 'deep'
ENDIF
IF keyword_set(johnson) then begin
    johnson_U = johnson_U[ii]
    johnson_B = johnson_B[ii]
    johnson_V = johnson_V[ii]
    johnson_R = johnson_R[ii]
    johnson_I = johnson_I[ii]
    mjohnson_U = mjohnson_U[ii]
    mjohnson_B = mjohnson_B[ii]
    mjohnson_V = mjohnson_V[ii]
    mjohnson_R = mjohnson_R[ii]
    mjohnson_I = mjohnson_I[ii]
    print, 'johnson'
ENDIF
IF keyword_set(ubv_deep) then begin
    ubv_deep_U = ubv_deep_U[ii]
    ubv_deep_B = ubv_deep_B[ii]
    ubv_deep_V = ubv_deep_V[ii]
    mubv_deep_U = mubv_deep_U[ii]
    mubv_deep_B = mubv_deep_B[ii]
    mubv_deep_V = mubv_deep_V[ii]
    print, 'uvbdeep'
ENDIF
IF keyword_set(alhambra) then begin
    filter1 = filter1[ii]
    filter2 = filter2[ii]
    filter3 = filter3[ii]
    filter4 = filter4[ii]
    filter5 = filter5[ii]
    filter6 = filter6[ii]
    filter7 = filter7[ii]
    filter8 = filter8[ii]
    filter9 = filter9[ii]
    filter10 = filter10[ii]
    filter11 = filter11[ii]
    filter12 = filter12[ii]
    filter13 = filter13[ii]
    filter14 = filter14[ii]
    filter15 = filter15[ii]
    filter16 = filter16[ii]
    filter17 = filter17[ii]
    filter18 = filter18[ii]
    filter19 = filter19[ii]
    filter20 = filter20[ii]
    mfilter1 = mfilter1[ii]
    mfilter2 = mfilter2[ii]
    mfilter3 = mfilter3[ii]
    mfilter4 = mfilter4[ii]
    mfilter5 = mfilter5[ii]
    mfilter6 = mfilter6[ii]
    mfilter7 = mfilter7[ii]
    mfilter8 = mfilter8[ii]
    mfilter9 = mfilter9[ii]
    mfilter10 = mfilter10[ii]
    mfilter11 = mfilter11[ii]
    mfilter12 = mfilter12[ii]
    mfilter13 = mfilter13[ii]
    mfilter14 = mfilter14[ii]
    mfilter15 = mfilter15[ii]
    mfilter16 = mfilter16[ii]
    mfilter17 = mfilter17[ii]
    mfilter18 = mfilter18[ii]
    mfilter19 = mfilter19[ii]
    mfilter20 = mfilter20[ii]
    print, 'alhambra'
ENDIF
IF keyword_set(flamex) then begin
    flamex_Bw = flamex_Bw[ii]
    flamex_R = flamex_R[ii]
    flamex_I = flamex_I[ii]
    flamex_J = flamex_J[ii]
    flamex_H = flamex_H[ii]
    flamex_Ks = flamex_Ks[ii]
    flamex_35 = flamex_35[ii]
    flamex_46 = flamex_46[ii]
    mflamex_Bw = mflamex_Bw[ii]
    mflamex_R = mflamex_R[ii]
    mflamex_I = mflamex_I[ii]
    mflamex_J = mflamex_J[ii]
    mflamex_H = mflamex_K[ii]
    mflamex_Ks = mflamex_Ks[ii]
    mflamex_35 = mflamex_35[ii]
    mflamex_46 = mflamex_46[ii]
    print, 'flamex'
ENDIF
IF keyword_set(cfhtls) then begin
    cfhtls_u = cfhtls_u[ii]
    cfhtls_g = cfhtls_g[ii]
    cfhtls_r = cfhtls_r[ii]
    cfhtls_i = cfhtls_i[ii]
    cfhtls_z = cfhtls_z[ii]
    mcfhtls_u = mcfhtls_u[ii]
    mcfhtls_g = mcfhtls_g[ii]
    mcfhtls_r = mcfhtls_r[ii]
    mcfhtls_i = mcfhtls_i[ii]
    mcfhtls_z = mcfhtls_z[ii]
    print, 'cfhtls'
ENDIF
IF keyword_set(hsc) then begin
    hsc_g = hsc_g[ii]
    hsc_r = hsc_r[ii]
    hsc_i = hsc_i[ii]
    hsc_z = hsc_z[ii]
    hsc_y = hsc_y[ii]
    mhsc_g = mhsc_g[ii]
    mhsc_r = mhsc_r[ii]
    mhsc_i = mhsc_i[ii]
    mhsc_z = mhsc_z[ii]
    mhsc_y = mhsc_y[ii]
    print, 'hsc'
ENDIF
IF keyword_set(lsst) then begin
    lsst_u = lsst_u[ii]
    lsst_g = lsst_g[ii]
    lsst_r = lsst_r[ii]
    lsst_i = lsst_i[ii]
    lsst_z = lsst_z[ii]
    lsst_y2 = lsst_y2[ii]
    lsst_y3 = lsst_y3[ii]
    lsst_y4 = lsst_y4[ii]
    mlsst_u = mlsst_u[ii]
    mlsst_g = mlsst_g[ii]
    mlsst_r = mlsst_r[ii]
    mlsst_i = mlsst_i[ii]
    mlsst_z = mlsst_z[ii]
    mlsst_y2 = mlsst_y2[ii]
    mlsst_y3 = mlsst_y3[ii]
    mlsst_y4 = mlsst_y4[ii]
    print, 'lsst'
ENDIF
IF keyword_set(euclid) then begin
    euclid_vis = euclid_vis[ii]
    euclid_y = euclid_y[ii]
    euclid_j = euclid_j[ii]
    euclid_h = euclid_h[ii]
    meuclid_vis = meuclid_vis[ii]
    meuclid_y = meuclid_y[ii]
    meuclid_j = meuclid_j[ii]
    meuclid_h = meuclid_h[ii]
    print, 'euclid'
ENDIF
IF keyword_set(irac) then begin
    irac_1 = irac_1[ii]
    irac_2 = irac_2[ii]
    irac_3 = irac_3[ii]
    irac_4 = irac_4[ii]
    mirac_1 = mirac_1[ii]
    mirac_2 = mirac_2[ii]
    mirac_3 = mirac_3[ii]
    mirac_4 = mirac_4[ii]
    print, 'irac'
ENDIF
IF keyword_set(wise) then begin
    wise_1 = wise_1[ii]
    wise_2 = wise_2[ii]
    wise_3 = wise_3[ii]
    wise_4 = wise_4[ii]
    mwise_1 = mwise_1[ii]
    mwise_2 = mwise_2[ii]
    mwise_3 = mwise_3[ii]
    mwise_4 = mwise_4[ii]
    print, 'wise'
ENDIF
IF keyword_set(wfirst) then begin
    wfirst_y = wfirst_y[ii]
    wfirst_j = wfirst_j[ii]
    wfirst_h = wfirst_h[ii]
    wfirst_k = wfirst_k[ii]
    mwfirst_y = mwfirst_y[ii]
    mwfirst_j = mwfirst_j[ii]
    mwfirst_h = mwfirst_h[ii]
    mwfirst_k = mwfirst_k[ii]
    print, 'wfirst'
ENDIF
IF keyword_set(twomass) then begin
    twomass_j = twomass_j[ii]
    twomass_h = twomass_h[ii]
    twomass_ks = twomass_ks[ii]
    mtwomass_j = mtwomass_j[ii]
    mtwomass_h = mtwomass_h[ii]
    mtwomass_ks = mtwomass_ks[ii]
    print, 'twomass'
ENDIF
IF keyword_set(candels) then begin
   candels_0 = candels_0[ii]
   candels_1 = candels_1[ii]
   candels_2 = candels_2[ii]
   candels_3 = candels_3[ii]
   candels_4 = candels_4[ii]
   candels_5 = candels_5[ii]
   candels_6 = candels_6[ii]
   candels_7 = candels_7[ii]
   candels_8 = candels_8[ii]
   candels_9 = candels_9[ii]
   mcandels_0 = mcandels_0[ii]
   mcandels_1 = mcandels_1[ii]
   mcandels_2 = mcandels_2[ii]
   mcandels_3 = mcandels_3[ii]
   mcandels_4 = mcandels_4[ii]
   mcandels_5 = mcandels_5[ii]
   mcandels_6 = mcandels_6[ii]
   mcandels_7 = mcandels_7[ii]
   mcandels_8 = mcandels_8[ii]
   mcandels_9 = mcandels_9[ii]
ENDIF
IF keyword_set(galz) then begin
   gz = gz[ii]
   zbin = zbin[ii]
   gzcentral = gzcentral[ii]
ENDIF


endif
   str.omag[0] = u
   str.omag[1] = g
   str.omag[2] = r
   str.omag[3] = i
   str.omag[4] = z
   str.amag[0] = mu
   str.amag[1] = mg
   str.amag[2] = mr
   str.amag[3] = mi
   str.amag[4] = mz
   tmp = k_sdss_bell(str.amag)
   print, 'k_sdss_bell'
   tmp1 = fltarr(ng)
   tmp1(*) = tmp(0,*)
   str.mstar = tmp1
   str.z = zz
   str.ra = ra
   str.dec = dec
   str.px = px
   str.py = py
   str.pz = pz
   str.vx = vx
   str.vy = vy
   str.vz = vz
   str.edge = 0
   str.siglos = siglos
   str.halopx = halopx
   str.halopy = halopy
   str.halopz = halopz
   str.halovx = halovx
   str.halovy = halovy
   str.halovz = halovz
   str.haloid = haloid
   str.m200 = m200
   str.r200 = r200
   str.rhalo = rhalo
   str.ngals = ngals
   str.haloz = haloz
   str.halora = halora
   str.halodec = halodec
   str.d8 = d8
   str.d8_drawn = td8
   str.central = central
   str.ecatid = ecatid
   str.nndist = nndist
   str.nnpercent = nnpercent
IF keyword_set(des) then begin
;    str.odesmag[0] = des_y
;    str.odesmag[1] = des_z 
;    str.odesmag[2] = des_g 
;    str.odesmag[3] = des_i 
;    str.odesmag[4] = des_r 
;    str.odesmag[5] = des_z2
;    str.adesmag[0] = mdes_y
;    str.adesmag[1] = mdes_z 
;    str.adesmag[2] = mdes_g 
;    str.adesmag[3] = mdes_i 
;    str.adesmag[4] = mdes_r 
;    str.adesmag[5] = mdes_z2
   str.odesmag[0] = des_u
   str.odesmag[1] = des_g
   str.odesmag[2] = des_r
   str.odesmag[3] = des_i
   str.odesmag[4] = des_z
   str.odesmag[5] = des_y
   str.adesmag[0] = mdes_u
   str.adesmag[1] = mdes_g
   str.adesmag[2] = mdes_r
   str.adesmag[3] = mdes_i
   str.adesmag[4] = mdes_z
   str.adesmag[5] = mdes_y
ENDIF
IF keyword_set(vista) then begin
    str.ovistamag[0] = wfcam_j
    str.ovistamag[1] = wfcam_h
    str.ovistamag[2] = wfcam_k
    str.ovistamag[3] = wfcam_y
    str.ovistamag[4] = wfcam_z
    str.avistamag[0] = mwfcam_j
    str.avistamag[1] = mwfcam_h
    str.avistamag[2] = mwfcam_k
    str.avistamag[3] = mwfcam_y
    str.avistamag[4] = mwfcam_z
    print, 'vista'
ENDIF
IF keyword_set(deep) then begin
    str.odeepmag[0] = deep_b
    str.odeepmag[1] = deep_r
    str.odeepmag[2] = deep_i
    str.adeepmag[0] = mdeep_b
    str.adeepmag[1] = mdeep_r
    str.adeepmag[2] = mdeep_i
    print, 'deep'
ENDIF
IF keyword_set(johnson) then begin
    str.ojohnsonmag[0] = johnson_U
    str.ojohnsonmag[1] = johnson_B
    str.ojohnsonmag[2] = johnson_V
    str.ojohnsonmag[3] = johnson_R
    str.ojohnsonmag[4] = johnson_I
    str.ajohnsonmag[0] = mjohnson_U
    str.ajohnsonmag[1] = mjohnson_B
    str.ajohnsonmag[2] = mjohnson_V
    str.ajohnsonmag[3] = mjohnson_R
    str.ajohnsonmag[4] = mjohnson_I
    print, 'johnson'
ENDIF
IF keyword_set(ubv_deep) then begin
    str.oubv_deepmag[0] = ubv_deep_U
    str.oubv_deepmag[1] = ubv_deep_B
    str.oubv_deepmag[2] = ubv_deep_V
    str.aubv_deepmag[0] = mubv_deep_U
    str.aubv_deepmag[1] = mubv_deep_B
    str.aubv_deepmag[2] = mubv_deep_V
    print, 'deep'
ENDIF
IF keyword_set(alhambra) then begin
    str.oalhambramag[0] = filter1
    str.oalhambramag[1] = filter2
    str.oalhambramag[2] = filter3
    str.oalhambramag[3] = filter4
    str.oalhambramag[4] = filter5
    str.oalhambramag[5] = filter6
    str.oalhambramag[6] = filter7
    str.oalhambramag[7] = filter8
    str.oalhambramag[8] = filter9
    str.oalhambramag[9] = filter10
    str.oalhambramag[10] = filter11
    str.oalhambramag[11] = filter12
    str.oalhambramag[12] = filter13
    str.oalhambramag[13] = filter14
    str.oalhambramag[14] = filter15
    str.oalhambramag[15] = filter16
    str.oalhambramag[16] = filter17
    str.oalhambramag[17] = filter18
    str.oalhambramag[18] = filter19
    str.oalhambramag[19] = filter20
    str.aalhambramag[0] = mfilter1
    str.aalhambramag[1] = mfilter2
    str.aalhambramag[2] = mfilter3
    str.aalhambramag[3] = mfilter4
    str.aalhambramag[4] = mfilter5
    str.aalhambramag[5] = mfilter6
    str.aalhambramag[6] = mfilter7
    str.aalhambramag[7] = mfilter8
    str.aalhambramag[8] = mfilter9
    str.aalhambramag[9] = mfilter10
    str.aalhambramag[10] = mfilter11
    str.aalhambramag[11] = mfilter12
    str.aalhambramag[12] = mfilter13
    str.aalhambramag[13] = mfilter14
    str.aalhambramag[14] = mfilter15
    str.aalhambramag[15] = mfilter16
    str.aalhambramag[16] = mfilter17
    str.aalhambramag[17] = mfilter18
    str.aalhambramag[18] = mfilter19
    str.aalhambramag[19] = mfilter20
    print, 'alhambra'
ENDIF
IF keyword_set(flamex) then begin
    str.oflamexmag[0] = flamex_Bw
    str.oflamexmag[1] = flamex_R
    str.oflamexmag[2] = flamex_I
    str.oflamexmag[3] = flamex_J
    str.oflamexmag[4] = flamex_H
    str.oflamexmag[5] = flamex_Ks
    str.oflamexmag[6] = flamex_36
    str.oflamexmag[7] = flamex_45
    str.aflamexmag[0] = mflamex_Bw
    str.aflamexmag[1] = mflamex_R
    str.aflamexmag[2] = mflamex_I
    str.aflamexmag[3] = mflamex_J
    str.aflamexmag[4] = mflamex_H
    str.aflamexmag[5] = mflamex_Ks
    str.aflamexmag[6] = mflamex_36
    str.aflamexmag[7] = mflamex_45
    print, 'flamex'
ENDIF
IF keyword_set(cfhtls) then begin
    str.ocfhtlsmag[0] = cfhtls_u
    str.ocfhtlsmag[1] = cfhtls_g
    str.ocfhtlsmag[2] = cfhtls_r
    str.ocfhtlsmag[3] = cfhtls_i
    str.ocfhtlsmag[4] = cfhtls_z
    str.acfhtlsmag[0] = mcfhtls_u
    str.acfhtlsmag[1] = mcfhtls_g
    str.acfhtlsmag[2] = mcfhtls_r
    str.acfhtlsmag[3] = mcfhtls_i
    str.acfhtlsmag[4] = mcfhtls_z
    print, 'cfhtls'
ENDIF
IF keyword_set(hsc) then begin
;  help, hsc_g, hsc_r, hsc_i, hsc_z, hsc_y
  str.ohscmag[0] = hsc_g
  str.ohscmag[1] = hsc_r
  str.ohscmag[2] = hsc_i
  str.ohscmag[3] = hsc_z
  str.ohscmag[4] = hsc_y
  str.ahscmag[0] = mhsc_g
  str.ahscmag[1] = mhsc_r
  str.ahscmag[2] = mhsc_i
  str.ahscmag[3] = mhsc_z
  str.ahscmag[4] = mhsc_y
  print, 'hsc'
ENDIF
IF keyword_set(lsst) then begin
  str.olsstmag[0] = lsst_u
  str.olsstmag[1] = lsst_g
  str.olsstmag[2] = lsst_r
  str.olsstmag[3] = lsst_i
  str.olsstmag[4] = lsst_z
  str.olsstmag[5] = lsst_y2
  str.olsstmag[6] = lsst_y3
  str.olsstmag[7] = lsst_y4
  str.alsstmag[0] = mlsst_u
  str.alsstmag[1] = mlsst_g
  str.alsstmag[2] = mlsst_r
  str.alsstmag[3] = mlsst_i
  str.alsstmag[4] = mlsst_z
  str.alsstmag[5] = mlsst_y2
  str.alsstmag[6] = mlsst_y3
  str.alsstmag[7] = mlsst_y4
  print, 'lsst'
ENDIF
IF keyword_set(euclid) then begin
  str.oeuclidmag[0] = euclid_vis
  str.oeuclidmag[1] = euclid_y
  str.oeuclidmag[2] = euclid_j
  str.oeuclidmag[3] = euclid_h
  str.aeuclidmag[0] = meuclid_vis
  str.aeuclidmag[1] = meuclid_y
  str.aeuclidmag[2] = meuclid_j
  str.aeuclidmag[3] = meuclid_h
  print, 'euclid'
ENDIF
IF keyword_set(irac) then begin
  str.oiracmag[0] = irac_1
  str.oiracmag[1] = irac_2
  str.oiracmag[2] = irac_3
  str.oiracmag[3] = irac_4
  str.airacmag[0] = mirac_1
  str.airacmag[1] = mirac_2
  str.airacmag[2] = mirac_3
  str.airacmag[3] = mirac_4
  print, 'irac'
ENDIF
IF keyword_set(wise) then begin
  str.owisemag[0] = wise_1
  str.owisemag[1] = wise_2
  str.owisemag[2] = wise_3
  str.owisemag[3] = wise_4
  str.awisemag[0] = mwise_1
  str.awisemag[1] = mwise_2
  str.awisemag[2] = mwise_3
  str.awisemag[3] = mwise_4
  print, 'wise'
ENDIF
IF keyword_set(wfirst) then begin
  str.owfirstmag[0] = wfirst_y
  str.owfirstmag[1] = wfirst_j
  str.owfirstmag[2] = wfirst_h
  str.owfirstmag[3] = wfirst_k
  str.awfirstmag[0] = mwfirst_y
  str.awfirstmag[1] = mwfirst_j
  str.awfirstmag[2] = mwfirst_h
  str.awfirstmag[3] = mwfirst_k
  print, 'wfirst'
ENDIF
IF keyword_set(twomass) then begin
  str.otwomassmag[0] = twomass_j
  str.otwomassmag[1] = twomass_h
  str.otwomassmag[2] = twomass_ks
  str.atwomassmag[0] = mtwomass_j
  str.atwomassmag[1] = mtwomass_h
  str.atwomassmag[2] = mtwomass_ks
  print, 'twomass'
ENDIF
IF keyword_set(candels) then begin
  str.ocandelsmag[0] = candels_0
  str.ocandelsmag[1] = candels_1
  str.ocandelsmag[2] = candels_2
  str.ocandelsmag[3] = candels_3
  str.ocandelsmag[4] = candels_4
  str.ocandelsmag[5] = candels_5
  str.ocandelsmag[6] = candels_6
  str.ocandelsmag[7] = candels_7
  str.ocandelsmag[8] = candels_8
  str.ocandelsmag[9] = candels_9
  str.acandelsmag[0] = mcandels_0
  str.acandelsmag[1] = mcandels_1
  str.acandelsmag[2] = mcandels_2
  str.acandelsmag[3] = mcandels_3
  str.acandelsmag[4] = mcandels_4
  str.acandelsmag[5] = mcandels_5
  str.acandelsmag[6] = mcandels_6
  str.acandelsmag[7] = mcandels_7
  str.acandelsmag[8] = mcandels_8
  str.acandelsmag[9] = mcandels_9
  print, 'candels'
ENDIF
IF keyword_set(galz) then begin
   str.gz = gz
   str.zbin = zbin
   str.gzcentral = gzcentral
ENDIF

end
