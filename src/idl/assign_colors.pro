pro assign_colors, path=path, matches=matches, declimit=declimit, $
                   thesecoeffs=thesecoeffs, omag=omag, amag=amag, $
                   des=des, vista=vista, deep=deep,johnson=johnson,$
                   ubv_deep=ubv_deep, dim=dim, alhambra=alhambra, $
                   flamex=flamex, cfhtls=cfhtls, hsc=hsc, $
		   lsst=lsst, euclid=euclid, wfirst=wfirst, twomass=twomass, $
		   irac=irac, wise=wise, candels=candels, cut=cut
if not keyword_set(path) then path ='/data/risa/sdss/galcats/'
print, path
print, 'assigning colors'
;rdfloat, Path+'gal_ginfo1.dat', Id, C1, C2, C3,  Msim, Ra, Dec, Zz,
;central
if (cut eq 0) then begin 
   if (dim eq 0) then rdfloat, Path+'gal_ginfo1.dat', Id, Msim, Ra, Dec, Zz, central
   if (dim eq 1) then rdfloat, Path+'gal_ginfo1_dim.dat', Id, Msim, Ra, Dec, Zz, central
endif else begin
   rdfloat, Path+'gal_ginfo1.dat.cut', Id, Msim, Ra, Dec, Zz, central
endelse


IF keyword_set(declimit) THEN BEGIN
    ii = where(dec LT declimit)
    help, ii
    id = id[ii]
    msim = msim[ii]
    ra = ra[ii]
    dec = dec[ii]
    zz = zz[ii]
    central = central[ii]
endif
Range, zz
magmin = min(msim)
;IF (magmin LT -18) THEN  sdss_gals = mrdfits('/data/working/smaller_dim_dr3.fit', 1) $
;ELSE sdss_gals = mrdfits('/data/working/smaller_dr3.fit', 1)

;sdss_gals = mrdfits('/data/working/coeffs.fit', 1)
;sdss_gals = mrdfits('~/projects/addgals/idl/coeffs.fit', 1)
sdss_gals = mrdfits('/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/cooper/combined_dr6_cooper.fit',1)
matches = sdss_gals[id-1]    
range, id
help, sdss_gals
IF keyword_set(thesecoeffs) THEN coeffs = thesecoeffs[*, id-1] $
ELSE coeffs = matches.coeffs
;coeffs = matches[id-1].coeffs  

outdat25=1
reconstruct_from_coeff, zz, coeffs, dat, morig=msim, change=change, outdat25=outdat25

;abcorr = [0.042,-0.036,-0.015,-0.013,0.002]
; switch from AB mags to SDSS bands
;for i=0,4 do dat.omag[i] = dat.omag[i]+abcorr[i]

check = 0
;temporary add to check
if check then begin
    read_miller_dat, rdat, /velim, /dr1
    add_absolute, rdat, coeffs=coeffs
    
    n = n_elements(msim)
    seed = 192
    newcoeffs = fltarr(3, n)
    for gi=0L,n-1 do begin
        ii = where(rdat.amag[2] gt Msim[gi]-0.1 and rdat.amag[2] lt Msim[gi]+0.1)
        if (ii[0] eq -1 and msim[gi] lt -22.3) then ii = where(rdat.amag[2] lt -22.3)
        if ii[0] eq -1 then ii = where(rdat.amag[2] gt Msim[gi]-0.2 and rdat.amag[2] lt Msim[gi]+0.2)
        if ii[0] eq -1 then ii = where(rdat.amag[2] gt Msim[gi]-0.3 and rdat.amag[2] lt Msim[gi]+0.3)
        
        nmatch = n_elements(ii)
;    print, gi, nmatch, ii[0]
        thisone = ii[floor(randomu(seed)*nmatch)]
        newcoeffs[*,gi] = coeffs[*, thisone]
    endfor
    
    reconstruct_from_coeff, zz, newcoeffs, dat, morig=msim
endif

for i=0,4 do begin
    ii = where(dat.omag gt 0 and dat.omag lt 100 and dat.amag[i] gt -26 and dat.amag[i] lt 0)
    help, dat, ii
    range, dat[ii].omag[i]
    range, dat[ii].amag[i]
endfor
if cut eq 0 then begin
   openw,outfile,path+'gal_ginfo.dat',/get_lun
endif else begin
   openw, outfile, path+'gal_ginfo.dat.cut',/get_lun
endelse
niceprintf, outfile, dat.omag(0), dat.omag(1), dat.omag(2), dat.omag(3), $
  dat.omag(4), dat.amag(0), dat.amag(1), dat.amag(2), dat.amag(3), $
  dat.amag(4), ra, dec, zz, central, id-1
close,/all

if cut eq 0 then begin
   openw,outfile,path+'gal_ginfo25.dat',/get_lun
endif else begin
   openw, outfile, path+'gal_ginfo25.dat.cut',/get_lun
endelse
niceprintf, outfile, outdat25.amag(0), outdat25.amag(1), outdat25.amag(2), $
  outdat25.amag(3), outdat25.amag(4)
close, /all

IF(keyword_set(des)) then begin
  reconstruct_from_coeff, zz, coeffs, DESdat, change=change, instrument='des'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_des.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_des.dat.cut',/get_lun
  endelse
  niceprintf, outfile, DESdat.omag(0), DESdat.omag(1), DESdat.omag(2), DESdat.omag(3), $
	DESdat.omag(4), DESdat.omag(5), DESdat.amag(0), DESdat.amag(1), DESdat.amag(2), $
	DESdat.amag(3),DESdat.amag(4), DESdat.amag(5)
  close,/all
  print, 'DES'
ENDIF
IF(keyword_set(vista)) then begin
  reconstruct_from_coeff, zz, coeffs, VISTAdat, change=change, instrument='vista'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_vista.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_vista.dat.cut',/get_lun
  endelse

  niceprintf, outfile, VISTAdat.omag(0), VISTAdat.omag(1), VISTAdat.omag(2), VISTAdat.omag(3), $
	VISTAdat.omag(4), VISTAdat.amag(0), VISTAdat.amag(1), VISTAdat.amag(2), $
	VISTAdat.amag(3), VISTAdat.amag(4)
  close,/all
  print, 'VISTA'
ENDIF
IF(keyword_set(deep)) then begin
  reconstruct_from_coeff, zz, coeffs, DEEPdat, change=change, instrument='deep'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_deep.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_deep.dat.cut',/get_lun
  endelse

  niceprintf, outfile, DEEPdat.omag(0), DEEPdat.omag(1), DEEPdat.omag(2), $
        DEEPdat.amag(0), DEEPdat.amag(1), DEEPdat.amag(2)
  close,/all
ENDIF
IF(keyword_set(johnson)) then begin
  reconstruct_from_coeff, zz, coeffs, Johnsondat, change=change, instrument='johnson'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_johnson.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_johnson.dat.cut',/get_lun
  endelse
  niceprintf, outfile, Johnsondat.omag(0), Johnsondat.omag(1), Johnsondat.omag(2), $
              Johnsondat.omag(3), Johnsondat.omag(4), $
              Johnsondat.amag(0), Johnsondat.amag(1), Johnsondat.amag(2), $
              Johnsondat.amag(3), Johnsondat.amag(4)
  close,/all
ENDIF
IF(keyword_set(ubv_deep)) then begin
  reconstruct_from_coeff, zz, coeffs, ubv_deepdat, change=change, instrument='ubv_deep'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_uvb_deep.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_uvb_deep.dat.cut',/get_lun
  endelse

  niceprintf, outfile, ubv_deepdat.omag(0), ubv_deepdat.omag(1), ubv_deepdat.omag(2), $
              ubv_deepdat.amag(0), ubv_deepdat.amag(1), ubv_deepdat.amag(2)
  close,/all
ENDIF
IF(keyword_set(alhambra)) then begin
  reconstruct_from_coeff, zz, coeffs, alhambradat, change=change, instrument='alhambra'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_alhambra.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_alhambra.dat.cut',/get_lun
  endelse

  niceprintf, outfile, ALHAMBRAdat.omag(0), ALHAMBRAdat.omag(1), ALHAMBRAdat.omag(2), $
              ALHAMBRAdat.omag(3), ALHAMBRAdat.omag(4), ALHAMBRAdat.omag(5), $
              ALHAMBRAdat.omag(6), ALHAMBRAdat.omag(7), ALHAMBRAdat.omag(8), $
              ALHAMBRAdat.omag(9), ALHAMBRAdat.omag(10), ALHAMBRAdat.omag(11), $
              ALHAMBRAdat.omag(12), ALHAMBRAdat.omag(13), ALHAMBRAdat.omag(14), $
              ALHAMBRAdat.omag(15), ALHAMBRAdat.omag(16), ALHAMBRAdat.omag(17), $
              ALHAMBRAdat.omag(18), ALHAMBRAdat.omag(19), $
              ALHAMBRAdat.amag(0), ALHAMBRAdat.amag(1), ALHAMBRAdat.amag(2), $
              ALHAMBRAdat.amag(3), ALHAMBRAdat.amag(4), ALHAMBRAdat.amag(5), $
              ALHAMBRAdat.amag(6), ALHAMBRAdat.amag(7), ALHAMBRAdat.amag(8), $
              ALHAMBRAdat.amag(9), ALHAMBRAdat.amag(10), ALHAMBRAdat.amag(11), $
              ALHAMBRAdat.amag(12), ALHAMBRAdat.amag(13), ALHAMBRAdat.amag(14), $
              ALHAMBRAdat.amag(15), ALHAMBRAdat.amag(16), ALHAMBRAdat.amag(17), $
              ALHAMBRAdat.amag(18), ALHAMBRAdat.amag(19)
  close,/all
ENDIF
IF(keyword_set(flamex)) then begin
  reconstruct_from_coeff, zz, coeffs, flamexdat, change=change, instrument='flamex'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_flamex.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_flamex.dat.cut',/get_lun
  endelse

  niceprintf, outfile, flamexdat.omag(0),flamexdat.omag(1),flamexdat.omag(2), $
              flamexdat.omag(3),flamexdat.omag(4),flamexdat.omag(5), $
              flamexdat.omag(6),flamexdat.omag(7), $
              flamexdat.amag(0), flamexdat.amag(1), flamexdat.amag(2), $
              flamexdat.amag(3), flamexdat.amag(4), flamexdat.amag(5), $
              flamexdat.amag(6), flamexdat.amag(7)
  close,/all
ENDIF
IF(keyword_set(wise)) then begin
  reconstruct_from_coeff, zz, coeffs, wisedat, change=change, instrument='wise'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_wise.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_wise.dat.cut',/get_lun
  endelse

  niceprintf, outfile, wisedat.omag(0),wisedat.omag(1),wisedat.omag(2), $
              wisedat.omag(3), wisedat.amag(0), wisedat.amag(1), $
	      wisedat.amag(2), wisedat.amag(3)
  close, /all
ENDIF
IF(keyword_set(cfhtls)) then begin
  reconstruct_from_coeff, zz, coeffs, cfhtlsdat, change=change, instrument='cfhtls'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_cfhtls.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_cfhtls.dat.cut',/get_lun
  endelse

  niceprintf, outfile, cfhtlsdat.omag(0),cfhtlsdat.omag(1),cfhtlsdat.omag(2), $
              cfhtlsdat.omag(3),cfhtlsdat.omag(4), $
              cfhtlsdat.amag(0), cfhtlsdat.amag(1), cfhtlsdat.amag(2), $
              cfhtlsdat.amag(3), cfhtlsdat.amag(4)
  close,/all
ENDIF
IF(keyword_set(hsc)) then begin
  reconstruct_from_coeff, zz, coeffs, hscdat, change=change, instrument='hsc'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_hsc.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_hsc.dat.cut',/get_lun
  endelse

  niceprintf, outfile, hscdat.omag(0),hscdat.omag(1),hscdat.omag(2), $
              hscdat.omag(3),hscdat.omag(4), $
              hscdat.amag(0), hscdat.amag(1), hscdat.amag(2), $
              hscdat.amag(3), hscdat.amag(4)
  close,/all
ENDIF
IF(keyword_set(lsst)) then begin
  reconstruct_from_coeff, zz, coeffs, lsstdat, change=change, instrument='lsst'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_lsst.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_lsst.dat.cut',/get_lun
  endelse

  niceprintf, outfile, lsstdat.omag(0),lsstdat.omag(1),lsstdat.omag(2), $
              lsstdat.omag(3),lsstdat.omag(4), lsstdat.omag(5),$
	      lsstdat.omag(6),lsstdat.omag(7),$
              lsstdat.amag(0), lsstdat.amag(1), lsstdat.amag(2), $
              lsstdat.amag(3), lsstdat.amag(4),lsstdat.amag(5), $
	      lsstdat.amag(6), lsstdat.amag(7)
  close,/all
ENDIF
IF(keyword_set(euclid)) then begin
  reconstruct_from_coeff, zz, coeffs, eucliddat, change=change, instrument='euclid'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_euclid.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_euclid.dat.cut',/get_lun
  endelse

  niceprintf, outfile, eucliddat.omag(0),eucliddat.omag(1),eucliddat.omag(2), $
              eucliddat.omag(3), $
              eucliddat.amag(0), eucliddat.amag(1), eucliddat.amag(2), $
              eucliddat.amag(3)
  close,/all
ENDIF
IF(keyword_set(irac)) then begin
  reconstruct_from_coeff, zz, coeffs, iracdat, change=change, instrument='irac'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_irac.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_irac.dat.cut',/get_lun
  endelse

  niceprintf, outfile, iracdat.omag(0),iracdat.omag(1),iracdat.omag(2), $
              iracdat.omag(3), $
              iracdat.amag(0), iracdat.amag(1), iracdat.amag(2), $
              iracdat.amag(3)
  close,/all
ENDIF
IF(keyword_set(wfirst)) then begin
  reconstruct_from_coeff, zz, coeffs, wfirstdat, change=change, instrument='wfirst'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_wfirst.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_wfirst.dat.cut',/get_lun
  endelse

  niceprintf, outfile, wfirstdat.omag(0),wfirstdat.omag(1),wfirstdat.omag(2), $
              wfirstdat.omag(3), $
              wfirstdat.amag(0), wfirstdat.amag(1), wfirstdat.amag(2), $
              wfirstdat.amag(3)
  close,/all
ENDIF
IF(keyword_set(twomass)) then begin
  reconstruct_from_coeff, zz, coeffs, twomassdat, change=change, instrument='twomass'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_twomass.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_twomass.dat.cut',/get_lun
  endelse

  niceprintf, outfile, twomassdat.omag(0),twomassdat.omag(1),twomassdat.omag(2), $
              twomassdat.amag(0), twomassdat.amag(1), twomassdat.amag(2)
  close,/all
ENDIF
IF(keyword_set(candels)) then begin
  print, 'Starting CANDELS'
  reconstruct_from_coeff, zz, coeffs, candelsDAT, change=change, instrument='candels'
  if cut eq 0 then begin
     openw,outfile,path+'gal_ginfo_candels.dat',/get_lun
  endif else begin
     openw, outfile, path+'gal_ginfo_candels.dat.cut',/get_lun
  endelse
  
  niceprintf, outfile,candelsDAT.omag(0), candelsDAT.omag(1),candelsDAT.omag(2), $
              candelsDAT.omag(3), candelsDAT.omag(4),candelsDAT.omag(5), $
              candelsDAT.omag(6), candelsDAT.omag(7),candelsDAT.omag(8), $
              candelsDAT.omag(9),$
              candelsDAT.amag(0), candelsDAT.amag(1),candelsDAT.amag(2), $
              candelsDAT.amag(3), candelsDAT.amag(4),candelsDAT.amag(5), $
              candelsDAT.amag(6), candelsDAT.amag(7),candelsDAT.amag(8), $
              candelsDAT.amag(9)
  close,/all
ENDIF


end
