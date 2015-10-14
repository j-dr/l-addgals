pro masstolight, g, halos, slow=slow, mw=mw, art=art
;jj = where(g.ngals ge 8)
;big = g[jj]
;lum20 = fltarr(n)
;this is in units of Lstar
;lumr = 10^(-0.4*(g.amag[2]+20.44))
;this is in units of Lsun
lumr = 10^(-0.4*(g.amag[2]-4.76))
s1=create_struct('haloid', 0L, 'm200', 0.0, 'ngals', 0L, 'r200', 0.0, $
                 'lumtot', 0.0, 'lum20', 0.0, 'lcent', 0.0, 'siglos', 0.0, $
                 'halopx', 0.0, 'halopy', 0.0, 'halopz', 0.0, $
                 'halovx', 0.0, 'halovy', 0.0, 'halovz', 0.0, $
                 'ra', 0.0, 'dec', 0.0, 'z', 0.0, 'n19', 0, 'n20', 0, 'n21', 0, $
                 'n22', 0,'edge', 0)

;g.haloid = g.hostid
if not keyword_set(mw) then begin
    print, 'assuming not mw'
    if keyword_set(art) then id = rem_dup(g.hostid) else id = rem_dup(g.haloid)
    n = n_elements(id)
    halos=replicate(s1,n)
    Print, n
    if keyword_set(art) then halos.haloid = g[id].hostid else     halos.haloid = g[id].haloid
;    PRINT, 'HALOID RANGE'
;    RANGE, HALOS.HALOID
    help, rem_dup(halos.haloid)
    help, rem_dup(g.haloid)
;    help, rem_dup(g.hostid)

    if keyword_set(art) then begin
        for i=0L,n-1 do begin
            inh = where(g.hostid eq halos[i].haloid)
           print, inh, where(g.haloid eq halos[i].haloid)
            ;sinh = where(g.hostid eq halos.haloid and g.central eq 0))
            halos[i].ngals = n_elements(inh)
            halos[i].lcent =  lumr(where(g.haloid eq halos[i].haloid))

            halos[i].lumtot =  lumr(inh)
            if(i/1000 eq 0) then print, i, inh[0], halos[i].ngals, halos[i].lcent/halos[i].lumtot
        endfor 
        return
    endif     else begin
        halos.ngals = g[id].ngals
    endelse
    print, 'got n'
    if not keyword_set(art) then begin
        print, 'assuming hv'
        halos.m200 = g[id].m200
        halos.r200 = g[id].r200
        halos.siglos = g[id].siglos
        halos.ra = g[id].halora
        halos.dec = g[id].halodec
        halos.z = g[id].haloz
        halos.halopx = g[id].halopx
        halos.halopy = g[id].halopy
        halos.halopz = g[id].halopz
        inh = where(g.rhalo lt g.r200, counth)
        in19 = where(g.rhalo lt g.r200 and g.amag[2] lt -19, count19)
        in20 = where(g.rhalo lt g.r200 and g.amag[2] lt -20, count20)
        in21 = where(g.rhalo lt g.r200 and g.amag[2] lt -21, count21)
        in22 = where(g.rhalo lt g.r200 and g.amag[2] lt -22, count22)
        print, 'in19:'
        print, in19
    endif else begin
        print, 'assuming art'
        halos.r200 = g[id].rhost
        halos.m200 = g[id].mhost
        inh = where(g.amag[2], counth)
        in19 = where(g.amag[2] lt -19, count19)
        in20 = where(g.amag[2] lt -20, count20)
        in21 = where(g.amag[2] lt -21, coutn21)
        in22 = where(g.amag[2] lt -22, count22)
    endelse
endif else begin
    add_tags, halos, ['m200', 'lumtot', 'r200', 'lum20'], ['0.0', '0.0', '0.0', '0.0'], halos2
    halos = halos2
    halos.m200 = halos.mvir
    halos.r200 = halos.rvir
    n = n_elements(halos)
    inh = where(g.lr  and g.amag[2] lt -19)
    in20 = where(g.lr and g.amag[2] lt -20)
endelse

if not keyword_set(slow) then begin
    if not keyword_set(art) then begin
;        if (inh[0] ge 0) then h = histogram(g[inh].haloid, min=0, max=max(g.haloid), rev=rev)
;        if (in19[0] ge 0) then h19 = histogram(g[in19].haloid, min=0, max=max(g.haloid), rev=rev19)
;        if (in20[0] ge 0) then h2 = histogram(g[in20].haloid, min=0, max=max(g.haloid), rev=rev2)
;        if (in21[0] ge 0) then h21 = histogram(g[in21].haloid, min=0, max=max(g.haloid), rev=rev21)
;        if (in22[0] ge 0) then h22 = histogram(g[in22].haloid, min=0, max=max(g.haloid), rev=rev22)
        if (counth gt 1) then h = histogram(g[inh].haloid, min=0, max=max(g.haloid), rev=rev)
        if (count19 gt 1) then h19 = histogram(g[in19].haloid, min=0, max=max(g.haloid), rev=rev19)
        if (count20 gt 1) then h2 = histogram(g[in20].haloid, min=0, max=max(g.haloid), rev=rev2)
        if (count21 gt 1) then h21 = histogram(g[in21].haloid, min=0, max=max(g.haloid), rev=rev21)
        if (count22 gt 1) then h22 = histogram(g[in22].haloid, min=0, max=max(g.haloid), rev=rev22)
;    endif else begin
        if (counth gt 1) then begin
        for i = 0L, n-1 do begin
            haloid = halos[i].haloid
            if rev[haloid] ne rev[haloid+1] then begin
                winh = rev[rev[haloid]:rev[haloid+1]-1]
;                halos[i].n19 = n_elements(win19)
                winh = inh[winh]
                halos[i].lumtot = total(lumr[winh]) ;
            endif
        endfor
        endif

        if (count19 gt 1) then begin
        for i = 0L, n-1 do begin
            haloid = halos[i].haloid
            if rev19[haloid] ne rev19[haloid+1] then begin
                win19 = rev19[rev19[haloid]:rev19[haloid+1]-1]
                halos[i].n19 = n_elements(win19)
            endif
        endfor 
        endif

        if (count20 gt 1) then begin
        for i = 0L, n-1 do begin
            haloid = halos[i].haloid
            if rev2[haloid] ne rev2[haloid+1] then begin
                win20 = rev2[rev2[haloid]:rev2[haloid+1]-1]
                                ;help, win20
                halos[i].n20 = n_elements(win20)
                win20 = in20[win20]
                halos[i].lum20 = total(lumr[win20]) ;       
                
            endif
;        print, i, halos[i].ngals, halos[i].lum20, halos[i].lum20/halos[i].lumtot
        endfor 
        endif

        if (count21 gt 1) then begin
        for i = 0L, n-1 do begin
            haloid = halos[i].haloid
            if rev21[haloid] ne rev21[haloid+1] then begin
                win21 = rev21[rev21[haloid]:rev21[haloid+1]-1]
                halos[i].n21 = n_elements(win21)
            endif
        endfor 
        endif

        if (count22 gt 1) then begin
        for i = 0L, n-1 do begin
            haloid = halos[i].haloid
            if rev22[haloid] ne rev22[haloid+1] then begin
                win22 = rev22[rev22[haloid]:rev22[haloid+1]-1]
                halos[i].n22 = n_elements(win22)
            endif
        endfor 
        endif

    endif
    if keyword_set(art) then begin
                                ;       h3 = histogram(g[in21].haloid, min=0, max=max(g.haloid), rev=rev3)
        
        h2 = histogram(g[in22].hostid, min=0, max=max(g.hostid), rev=rev2)
        h3 = histogram(g[in21].hostid, min=0, max=max(g.hostid), rev=rev3)
        h4 = histogram(g[in20].hostid, min=0, max=max(g.hostid), rev=rev4)
        for i = 0L, n-1 do begin
            haloid = halos[i].haloid
            if rev2[haloid] ne rev2[haloid+1] then begin
                win22 = rev2[rev2[haloid]:rev2[haloid+1]-1]
                halos[i].n22 = n_elements(win22)
            endif
            if rev3[haloid] ne rev3[haloid+1] then begin
                win21 = rev3[rev3[haloid]:rev3[haloid+1]-1]
                halos[i].n21 = n_elements(win21)
            endif
            if rev4[haloid] ne rev4[haloid+1] then begin
                win20 = rev4[rev4[haloid]:rev4[haloid+1]-1]
                halos[i].n20 = n_elements(win20)
            endif
        endfor
    endif
endif

;endif else begin
;    for i = 0L, n-1 do begin
;        tall = where(g.haloid eq halos[i].haloid and g.rhalo lt g.r200)
;        t20 = where(g.haloid eq halos[i].haloid and g.rhalo lt g.r200 and g.amag[2] lt -20)
;        t21 = where(g.haloid eq halos[i].haloid and g.rhalo lt g.r200 and g.amag[2] lt -21)
;        t22 = where(g.haloid eq halos[i].haloid and g.rhalo lt g.r200 and g.amag[2] lt -22)
 ;       IF KEYWORD_SET(ART) THEN BEGIN
;            tall = where(g.haloid eq halos[i].haloid and g.rhalo lt g.r200)
;            t20 = where(g.haloid eq halos[i].haloid and g.rhalo lt g.r200 and g.amag[2] lt -20)
 ;           t21 = where(g.haloid eq halos[i].haloid and g.rhalo lt g.r200 and g.amag[2] lt -21)
 ;;           t22 = where(g.haloid eq halos[i].haloid and g.rhalo lt g.r200 and g.amag[2] lt -22)
 ;       ENDIF;;

;        if(tall[0] gt -1) then begin
;            halos[i].lumtot = total(lumr[tall])
;            halos[i].ngals = n_elements(tall)
;            halos[i].n20 = n_elements(t20)
;            halos[i].n21 = n_elements(t21)
 ;           halos[i].n22 = n_elements(t21)
  ;          print, halos[i].m200, halos[i].ngals
  ;      endif
  ;      if(t20[0] gt -1) then begin
  ;          halos[i].lum20 = total(lumr[t20])
; ;           halos[i].ngals = n_elements(t20)
 ;       endif
 ;       print, i, halos[i].ngals, n_elements(t20), halos[i].lum20, halos[i].lum20/halos[i].lumtot
;    endfor
;endelse

;ii = where(halos.lumtot gt 0)
;halos=halos[ii]

end
