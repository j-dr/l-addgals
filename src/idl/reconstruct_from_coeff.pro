pro reconstruct_from_coeff, zz, coeffs, outdat, morig=morig, instrument=instrument, $
	outdat25=outdat25, change=change
;	des=des, vista=vista,$
;	deep=deep, johnson=johnson, ubv_deep=ubv_deep, alhambra=alhambra, $
;	flamex=flamex, cfhtls=cfhtls, hsc=hsc, lsst=lsst, euclid=euclid,$
;	wfirst=wfirst

if not KEYWORD_SET(instrument) then instrument = 'SDSS'

;;;we need to apply an AB correction for SDSS, but everything else is fine, 
;;;so this is just a dummy array that is set to zero for everything but SDSS.  
abcorr = fltarr(50)

case (strupcase(instrument)) of
  'SDSS':begin
	filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', 'sdss_z0.par']
	z0 = 0.1
	; Reconstructed are like true AB mags.  
	; put them back into the SDSS bands
	abcorr=[-0.036, 0.012, 0.010, 0.028, 0.040]
  end
  'DES':begin
	filterlist=['des_u.par', 'des_g.par', 'des_r.par', 'des_i.par', 'des_z.par', 'des_Y.par']
;	filterlist=['DES_asahi_g.par', 'DES_asahi_r.par', 'DES_asahi_i.par', 'DES_asahi_z.par', 'DES_asahi_Y.par']
  	z0 = 0.1
  end

  'VISTA':begin
  	filterlist=['wfcam_Z.par', 'wfcam_Y.par', 'wfcam_J.par', 'wfcam_H.par', 'wfcam_K.par']
  	z0 = 0.1
  end
  'DEEP':begin
  	filterlist=['deep_B.par', 'deep_R.par', 'deep_I.par']
  	z0 = 0.0
  end
  'JOHNSON':begin
  	filterlist=['bessell_U.par', 'bessell_B.par', 'bessell_V.par', $
		    'bessell_R.par', 'bessell_I.par']
  	z0 = 0.0
  end
  'UBV_DEEP':begin
  	filterlist=['johnson_U.par', 'johnson_B.par', 'johnson_V.par']
  	z0 = 0.0
  end
  'ALHAMBRA':begin
  	filterlist=['alhambra1.par', 'alhambra2.par', 'alhambra3.par', 'alhambra4.par', 'alhambra5.par', $
              'alhambra6.par', 'alhambra7.par', 'alhambra8.par', 'alhambra9.par', 'alhambra10.par', $
              'alhambra11.par', 'alhambra12.par', 'alhambra13.par', 'alhambra14.par', 'alhambra15.par', $
              'alhambra16.par', 'alhambra17.par', 'alhambra18.par', 'alhambra19.par', 'alhambra20.par']
  	z0 = 0.0
  end
  'FLAMEX':begin
  	filterlist=['flamex_Bw.par', 'flamex_R.par', 'flamex_I.par', 'flamex_J.par', 'flamex_H.par', 'flamex_Ks.par', 'flamex_36.par', 'flamex_45.par']
  	z0 = 0.0
  end
  'CFHTLS':begin
	filterlist=['CFHTLS_u.par','CFHTLS_g.par','CFHTLS_r.par','CFHTLS_i.par','CFHTLS_z.par']
	z0 = 0.0
  end
  'HSC':begin
	filterlist=['hsc_g.par', 'hsc_r.par', 'hsc_i.par', 'hsc_z.par', 'hsc_y.par']
	z0 = 0.0
  end
  'LSST':begin
	filterlist=['lsst_u.par','lsst_g.par','lsst_r.par','lsst_i.par','lsst_z.par',$
		    'lsst_y2.par','lsst_y3.par','lsst_y4.par']
	z0 = 0.0
  end
  'EUCLID':begin
	filterlist=['euclid_vis.par', 'euclid_Y.par', 'euclid_J.par', 'euclid_H.par']
	z0 = 0.0
  end
  'IRAC':begin
        filterlist=['spitzer_irac_ch1.par', 'spitzer_irac_ch2.par', 'spitzer_irac_ch3.par', 'spitzer_irac_ch4.par']
        z0 = 0.0
  end
  'WISE':begin
        ;filterlist=['spitzer_wise_ch1.par', 'spitzer_wise_ch2.par', 'spitzer_wise_ch3.par', 'spitzer_wise_ch4.par']
        filterlist=['RSR-W1.dat', 'RSR-W2.dat','RSR-W3.dat','RSR-W4.dat']
        z0 = 0.0
  end
  'WFIRST':begin
	filterlist=['wfirst_Y.par', 'wfirst_J.par', 'wfirst_H.par', 'wfirst_K.par']
	z0 = 0.0
  end
  'TWOMASS':begin
        filterlist=['twomass_J.par', 'twomass_H.par', 'twomass_Ks.par']
	z0 = 0.0
  end
  'CANDELS':begin
        filterlist=[ 'CANDELS/ACS/f435w.WFC1.par', 'CANDELS/ACS/f606w.WFC1.par',$
                     'CANDELS/ACS/f775w.WFC1.par', 'CANDELS/ACS/f814w.WFC1.par',$
                     'CANDELS/ACS/f850lp.WFC1.par', 'CANDELS/WFC3/f275w.UVIS1.par',$
                     'CANDELS/WFC3/f336w.UVIS1.par', 'CANDELS/WFC3/f105w.IR.par',$
                     'CANDELS/WFC3/f125w.IR.par', 'CANDELS/WFC3/f160w.IR.par']
        z0 = 0.0
        print, 'CANDELS filter list set'
  end
  else:begin
	print, 'Unknown instrument name: '+instrument
	print, 'Not generating magnitudes.'
	return
  end
endcase
nbands = N_ELEMENTS(filterlist)

;;;pre-load the templates to speed up the running of k-correct later
;;;on
print, 'loading vmatrix'
k_load_vmatrix, vmatrix, lambda

;;;the structure for storing the output magnitudes
nobj= n_elements(zz)
s2=create_struct('omag',fltarr(nbands), 'amag', fltarr(nbands), 'z', 0.0)
outdat=replicate(s2,nobj)

;;;do the reconstruction in the observed frame
print, 'reconstructing maggies'
k_reconstruct_maggies,coeffs,zz,reconstruct_maggies0,vmatrix=vmatrix, lambda=lambda, filterlist=filterlist
range, reconstruct_maggies0

;In the rest frame
k_reconstruct_maggies,coeffs,replicate(z0, nobj),reconstruct_maggies,vmatrix=vmatrix, lambda=lambda, filterlist=filterlist

; get kcorrection
kcorrect=reconstruct_maggies/reconstruct_maggies0
kcorrect=2.5*alog10(kcorrect)

;;;get the distance modulus
dm = lf_distmod(zz)

;;;go from maggies to magnitudes
outdat.omag= 22.5-2.5*alog10(reconstruct_maggies0/1.e-9)
for i=0,nbands-1 do outdat.amag[i]=outdat.omag[i]-dm-kcorrect[i, *]


IF KEYWORD_SET(outdat25) then begin
  ;;;In Sarah's 0.25 band
  k_reconstruct_maggies,coeffs,replicate(0.25, nobj),reconstruct_maggies,vmatrix=vmatrix, lambda=lambda, filterlist=filterlist
  kcorrect=reconstruct_maggies/reconstruct_maggies0
  kcorrect=2.5*alog10(kcorrect)
  ss1 = create_struct('amag', fltarr(nbands))
  outdat25 = replicate(ss1, nobj)
  for i = 0, nbands - 1 do begin
    outdat25.amag[i] = outdat.omag[i]-dm-kcorrect[i,*]
  endfor
ENDIF

for i=0,nbands-1 do outdat.omag(i)=outdat.omag(i)-abcorr(i)

outdat.z = zz

;;;we have a slight perturbation to get the magnitudes back to the original SDSS mags.
if (keyword_set(morig) or keyword_set(change)) then begin
    if not keyword_set(change) then change = morig-outdat.amag[2]
    for i=0,nbands-1 do begin
        outdat.omag[i] = outdat.omag[i]+change
        outdat.amag[i] = outdat.amag[i]+change
        if KEYWORD_SET(outdat25) then outdat25.amag[i] += change
    endfor
endif

end

