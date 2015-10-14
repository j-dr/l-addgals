function trim_gals, g, h

nobj = N_ELEMENTS(g)
;;; This is our new, trimed galaxy structure
s1=create_struct('id',0L,'ecatid', 0L, 'omag',fltarr(5),$
                 'amag',fltarr(5),'z',0.0,'ra',0.0,'dec',0.0,$
                 'px',0.0,'py',0.0,'pz',0.0,'vx',0.0,'vy',0.0,'vz',0.0,$
                 'edge', 0, 'haloid',0L,'m200',0.0,'ngals',0,'r200',0.0,$
                 'rhalo',0.0,'central',0,'mstar',0., 'd8', 0., 'nndist', 0., $
                 'nnpercent', 0.)
gals=replicate(s1,nobj)
gals(*).id = g(*).id
gals(*).ecatid = g(*).ecatid
gals(*).omag(*) = g(*).omag(*)
gals(*).amag(*) = g(*).amag(*)
gals(*).z = g(*).z
gals(*).ra = g(*).ra
gals(*).dec = g(*).dec
gals(*).px = g(*).px
gals(*).py = g(*).py
gals(*).pz = g(*).pz
gals(*).vx = g(*).vx
gals(*).vy = g(*).vy
gals(*).vz = g(*).vz
gals(*).edge = 0
gals(*).haloid = g(*).haloid
gals(*).m200 = g(*).m200
gals(*).ngals = g(*).ngals
gals(*).r200 = g(*).r200
gals(*).rhalo = g(*).rhalo
gals(*).central = g(*).central
gals(*).mstar = g(*).mstar
gals(*).d8 = g(*).d8
gals(*).nndist = g(*).nndist
gals(*).nnpercent = g(*).nnpercent

;;;--- a lookup table to quickly find halos based on haloid
HaloIdTab = lonarr(max(h.haloid)+1)
NHalos = N_ELEMENTS(h)
for i = 0L, NHalos - 1 do begin
  HaloIdTab(h(i).haloid) = i
endfor

;;;--- loop through galaxies, removing gals with separation larger
;;;    than r200 from the ngals count
h(*).ngals = 0
for i = 0L, nobj-1 do begin
  if (gals(i).rhalo ge gals(i).r200) then begin
    gals(i).ngals = -99
  endif else begin
    ThisHalo = HaloIdTab(gals(i).haloid)
    h(ThisHalo).ngals++
    if (gals(i).central) then h(ThisHalo).lcent = gals(i).amag(2)
  endelse
endfor

;;;--- transfer the new ngals information to the galaxies
for i = 0L, nobj-1 do begin
  if (gals(i).rhalo lt gals(i).r200) then begin
    ThisHalo = HaloIdTab(gals(i).haloid)
    gals(i).ngals = h(ThisHalo).ngals
  endif
endfor

return, gals

END
