pro add_errors_to_sdss_simulation,gal,galn,zp=zp

    if not keyword_set(zp) then zp=[24.,27.,27.,26.,24.]

    ng=n_elements(gal)

    galnone=create_struct('id',0L,'tmag',fltarr(5),'omag',fltarr(5),$
       'omagerr',fltarr(5),$
       'ra',0.0,'dec',0.0,'z',0.0,'photoz',0.0,'haloid',0L,'rhalo',0.0,$
       'r200',0.0,'central',0)
    galn=replicate(galnone,ng)
    galn.id=gal.id
    galn.ra=gal.ra
    galn.dec=gal.dec
    galn.z=gal.z
    galn.tmag=gal.omag
    galn.haloid=gal.haloid
    galn.rhalo=gal.rhalo
    galn.r200=gal.r200
    galn.central=gal.central

    for i=0,4 do begin
      flux=10^((zp(i)-gal.omag(i))/2.5)
      err=sqrt(flux)
      flux_err=flux + randomn(seed,ng)*err
      galn.omag(i)=zp(i)-2.5*alog10(flux_err)
;      galn.omagerr(i)=2.5*alog10(1.0+err/flux) ;tim's original error
      galn.omagerr(i)=2.5*alog10(1.0+err*(1+randomn(seed,ng)*0.12)/flux)
    endfor

    zerr=0.03*randomn(seed,ng)
    galn.photoz=galn.z+zerr*(1.0+galn.z)

    j=where(galn.omagerr(3) le 0.1)
    galn=galn(j)

    return
    end
