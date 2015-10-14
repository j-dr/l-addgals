pro make_rcs_catalog, g, mag, hgn, lmag=lmag


  sky_mag = [22.9138, 22.2964, 21.9881, 21.5390] ;; my calibration to fit Eli's data
  sky_mag -= 0.1
  zps = [25.11, 24.80, 24.36, 22.83] ;; from sdss fluxcal webpage
  npixels = 11.2689
  flux_sky = 10.^(-0.4*(sky_mag-zps))
  nk = flux_sky
  exposure_time = [397.458, 555.225, 491.012, 195.870]
  exposure_time *= 1.1
  intrinsic_noise = [0.000936655, 0.00116858, 0.00116377, 9.38784e-05]

  ng=n_elements(g)

     if KEYWORD_SET(lmag) then begin
       hgnone=create_struct('id',0L,'tmag',fltarr(4),'lmag',fltarr(4),'omag',fltarr(4),$
                            'omagerr',fltarr(4),'amag', fltarr(4), $
                            'flux', fltarr(4), 'ivar', fltarr(4), $
                            'ra',0.0,'dec',0.0,'z',0.0,'photoz_Gaussian',0.0,'haloid',0L,'rhalo',0.0,$
                            'm200', 0.0, 'r200',0.0,'central',0, 'tra', 0., 'tdec', 0., $
                            'gamma1', 0., 'gamma2', 0., 'kappa', 0., 'mu', 0., 'size', 0., $
                            'tsize', 0., 'epsilon', fltarr(2), 'te', fltarr(2))
     endif else begin
       hgnone=create_struct('id',0L,'tmag',fltarr(4),'omag',fltarr(4),$
                            'omagerr',fltarr(4),'amag', fltarr(4), $
                            'flux', fltarr(4), 'ivar', fltarr(4), $
                            'ra',0.0,'dec',0.0,'z',0.0,'photoz_Gaussian',0.0,'haloid',0L,'rhalo',0.0,$
                            'm200', 0.0, 'r200',0.0,'central',0, 'tra', 0., 'tdec', 0.)
     endelse

     hgn=replicate(hgnone,ng)
     hgn.id=g.id
     hgn.ra=g.ra
     hgn.dec=g.dec
     hgn.z=g.z
     hgn.haloid=g.haloid
     hgn.rhalo=g.rhalo
     hgn.r200=g.r200
     hgn.m200=g.m200
     hgn.central=g.central
     hgn.tmag(0) = mag.omag(1)
     hgn.tmag(1) = mag.omag(2)
     hgn.tmag(2) = mag.omag(3)
     hgn.tmag(3) = mag.omag(4)
     hgn.amag(0) = mag.amag(1)
     hgn.amag(1) = mag.amag(2)
     hgn.amag(2) = mag.amag(3)
     hgn.amag(3) = mag.amag(4)

     if KEYWORD_SET(lmag) then begin
       hgn.lmag(0) = g.lmag(1)
       hgn.lmag(1) = g.lmag(2)
       hgn.lmag(2) = g.lmag(3)
       hgn.lmag(3) = g.lmag(4)
       hgn.gamma1 = g.gamma1
       hgn.gamma2 = g.gamma1
       hgn.kappa = g.kappa
       hgn.mu = g.mu
       hgn.size = g.size
       hgn.tsize = g.tsize
       hgn.epsilon = g.epsilon
       hgn.te = g.te
     endif


  for i = 0, 3 do begin
     if KEYWORD_SET(lmag) then $
        flux_gal = 10.^(-0.4*(hgn.lmag(i) - zps(i)))*exposure_time[i] else $
           flux_gal = 10.^(-0.4*(hgn.tmag(i) - zps(i)))*exposure_time[i]
     noise = sqrt(nk[i]*exposure_time[i]*npixels + flux_gal) + intrinsic_noise[i]*flux_gal
     flux_err = flux_gal + randomn(seed,ng)*noise
     if tag_exist(hgn, 'flux') then hgn.flux[i] = flux_err
     if tag_exist(hgn, 'ivar') then hgn.ivar[i] = 1./noise^2
     ibad = where(flux_err lt 0, count, comp=igood)
     hgn[igood].omag(i) = zps[i] - 2.5*alog10(flux_err[igood]/exposure_time[i])
     if (count gt 0) then hgn[ibad].omag(i) = 99
     hgn.omagerr(i) = 2.5*alog10(1.0 + noise/flux_gal)
  endfor

  zerr=0.03*randomn(seed,ng)
  hgn.photoz_Gaussian=hgn.z+zerr*(1.0+hgn.z)

;;;; Q: which of the below definitions of j is better??????
;    j=where(hgn.omagerr(1) le 0.1 and hgn.omagerr(2) le 0.1 and $
;            hgn.omagerr(3) le 0.1)
;    j=where(hgn.omagerr(3) le 0.1)
;    hgn=hgn(j)

    return

    end



