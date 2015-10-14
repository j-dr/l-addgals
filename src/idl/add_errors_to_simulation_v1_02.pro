pro add_errors_to_simulation_v1_02,hg,des,hgn,zp=zp

    if not keyword_set(zp) then zp=29.0

    ng=n_elements(hg)

    hgnone=create_struct('id',0L,'tmag',fltarr(6),'omag',fltarr(6),$
       'omagerr',fltarr(6),$
       'ra',0.0,'dec',0.0,'z',0.0,'photoz',0.0,'haloid',0L,'rhalo',0.0,$
       'r200',0.0,'central',0)
    hgn=replicate(hgnone,ng)
    hgn.id=hg.id
    hgn.ra=hg.ra
    hgn.dec=hg.dec
    hgn.z=hg.z
    hgn.tmag=des.omag
    hgn.haloid=hg.haloid
    hgn.rhalo=hg.rhalo
    hgn.r200=hg.r200
    hgn.central=hg.central

    for i=0,5 do begin
      flux=10^((zp-des.omag(i))/2.5)
      err=sqrt(flux)
      flux_err=flux + randomn(seed,ng)*err
      hgn.omag(i)=zp-2.5*alog10(flux_err)
      hgn.omagerr(i)=2.5*alog10(1.0+err/flux) 
    endfor

    zerr=0.03*randomn(seed,ng)
    hgn.photoz=hgn.z+zerr*(1.0+hgn.z)

;;;; Q: which of the below definitions of j is better??????
;    j=where(hgn.omagerr(1) le 0.1 and hgn.omagerr(2) le 0.1 and $
;            hgn.omagerr(3) le 0.1)
    j=where(hgn.omagerr(3) le 0.1)
    hgn=hgn(j)

    return

    end



