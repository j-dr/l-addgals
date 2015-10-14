
; the magnitdues we're going to look at
mag = ['18', '19', '20', '21']
type = ['all']
carr = [!blue, !green, !orange, !red]

; the simulation files
simbase = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla-1/Catalog_v2.0/wp_cat/wp_galaxies_wp_cuts.fit_'
obsbase = '/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr7/wp/wp_sdss_dr7safe0_'

; loop over the types
for i = 0, N_ELEMENTS(type) - 1 do begin

  ; loop over the magnitudes
  for j = 0L, N_ELEMENTS(mag) - 1 do begin

    obsfile = obsbase+mag[j]+'_'+type[i]+'.dat'
    simfile = simbase+mag[j]+'_'+type[i]+'.dat'

    ; read the data
    rdfloat, obsfile, ologrp, orpavg, owp, onpair
    orp = 10.^ologrp
    rdfloat, simfile, slogrp, srpavg, swp, snpair
    srp = 10.^slogrp

    ; setup an empty plot
    if (j eq 0) then begin
      plot, orp, owp, /xlog, /ylog, yrange = [0.1,1e4], $
	xtitle = 'r!Dp', ytitle = 'w!Dp!N(r!Dp!N)', /nodat
    endif

    ; plot the data
    oplot, orp, owp, color = carr[j]
    oplot, srp, swp, color = carr[j], linestyle = 2

  endfor
  legend, 'Mr < -'+mag, carr, linestyle=0, /right
endfor

end
