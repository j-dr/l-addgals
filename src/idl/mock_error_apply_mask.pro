pro mock_error_apply_mask,tmag,mask,oflux,ofluxerr,omag,omagerr,seed=seed, point_source=point_source

lnscat = [0.1, 0.1, 0.1, 0.1, 0.1]
mask.lim_exptime[4] = 21.459
mask.lim_limmag[4] = 728.00

if n_params() eq 0 then begin
    print,'syntax- mock_error_apply_mask,tmag,mask,oflux,ofluxerr,omag,omagerr,seed=seed'
    print,' takes an array of tmags and mask structure with lim_exptime and lim_limmag, returns arrays of omags, omagerrs'
    return
endif

if n_elements(seed) eq 0 then seed=systime(/seconds)

nmag = N_ELEMENTS(tmag[*,0])

if (nmag gt n_elements(mask[0].lim_exptime)) then begin
    print,'Error: not enough magnitude errors sepcified'
    return
endif

ngal=n_elements(tmag[0,*])

zp=22.5

oflux=fltarr(nmag,ngal)
ofluxerr=oflux
omag=oflux
omagerr=oflux

offset = 0.0
if KEYWORD_SET(point_source) then begin
    offset = 0.6
    lnscat /= 4.0
endif

flux1_lim = fltarr(ngal)
fsky1 = fltarr(ngal)
tflux = fltarr(ngal)
flux = fltarr(ngal)
noise = fltarr(ngal)
for i=0l,nmag-1 do begin
  ii = where(mask.lim_exptime[0] gt 0.0, count,comp=jj) ;;cut objects not observed in this band

  ;;;calculate 1-second flux of sky
  flux1_lim[ii] = 10.^((mask[ii].lim_exptime[i]-zp)/(-2.5)) > 120/mask[ii].lim_exptime[i] ;; don't go negative...
  fsky1[ii] = (flux1_lim[ii]^2.*mask[ii].lim_exptime[i])/100. - flux1_lim[ii]

    tflux[ii] = mask[ii].lim_exptime[i] * 10.^((reform(tmag[i,ii])-offset-zp)/(-2.5))

    noise[ii] = exp(alog(sqrt(fsky1[ii]*mask[ii].lim_exptime[i] + tflux[ii])) + lnscat[i]*randomn(seed,count))

    flux[ii] = tflux[ii] + noise[ii]*randomn(seed,count)

    oflux[i,ii] = flux[ii] / mask[ii].lim_exptime[i]
    ofluxerr[i,ii] = noise[ii]/mask[ii].lim_exptime[i]
    if (count lt ngal) then begin
      oflux[i,jj] = 0.0
      ofluxerr[i,jj] = 0.0
    endif

    omag[i,*] = 22.5-2.5*alog10(oflux[i,*])
    omagerr[i,*] = (2.5/alog(10.))*(ofluxerr[i,*]/oflux[i,*])

    bad=where(finite(omag[i,*]) eq 0,nbad)
    if (nbad gt 0) then begin
        omag[i,bad] = 99.0
        omagerr[i,bad] = 99.0
    endif

endfor

return
end

