pro add_to_hod_nonperiodic, h, g, add_centrals=add_centrals, omegam=omegam, omegal=omegal, min_mass = min_mass, mag_ind = mag_ind, mass_def=mass_def

IF NOT KEYWORD_SET(omegam) then omegam = 0.25
IF NOT KEYWORD_SET(omegal) then omegal = 1.0 - omegam
IF NOT KEYWORD_SET(min_mass) then min_mass = 0
IF NOT KEYWORD_SET(mag_ind) then mag_ind = 2

case (strupcase(mass_def)) of
  'M200':begin
    mass = h.m200
    rad = h.r200
  end
  'M200b':begin
    mass = h.m200b
    rad = h.r200b
  end
  'M500':begin
    mass = h.m500
    rad = h.r500
  end
  'MVIR':begin
    mass = h.mvir
    rad = h.rvir
  end
  ELSE:begin
    print, "No mass definition specified.  Assuming m200c."
    mass = h.m200
    rad = h.r200
  end
endcase

print, 'mass definition used: '+mass_def

Q0 = 2.0
QZ0 = 0.1
Q1 = -1.0

Mr = g.amag(mag_ind) + Q0*(1 + Q1*(g.z - QZ0))*(g.z - QZ0)
;Mr = g.amag(mag_ind) + 1.3*(g.z - 0.1)

g.haloid = -1
if (TAG_EXIST(g, 'm200')) then g.m200 = 0.
g.r200 = 0.
if (TAG_EXIST(g, 'ngals')) then g.ngals = 0
if KEYWORD_SET(add_centrals) then g.central = 0
pid = lindgen(N_ELEMENTS(h))
if (TAG_EXIST(h, 'pid')) then cc = where(h.pid lt 0)
if (TAG_EXIST(h, 'host_haloid')) then cc = where(h.host_haloid eq h.haloid)

ng = N_ELEMENTS(g)
print, 'number of gals: ', ng
print, 'finding closest halos'
;;;find the closest halos
hcx = h[cc].halopx
hcy = h[cc].halopy
hcz = h[cc].halopz
hid=MATCH_3D(g.px,g.py,g.pz,hcx,hcy,hcz, 15.0, match_distance=r_closest)
print, 'size(hid):', size(hid)
good = where(hid ge 0)
print, 'size(good):', size(good)
hid[good] = cc[hid[good]]


print, 'looping through galaxies and adding info'
;;;loop through the galaxies to add the informaiton
for i = 0L, ng - 1 do begin
  if (r_closest[i] gt 5) then hid[i] = -1
  if (hid[i] eq -1) then begin
    g[i].rhalo = 9999.
    continue
  endif
  g[i].haloid = h[hid[i]].haloid
  g[i].rhalo = r_closest[i]
  if (g[i].rhalo lt rad[hid[i]]) then begin
     if (TAG_EXIST(g, 'm200')) then g[i].m200 = mass[hid[i]]
     if (TAG_EXIST(g, 'r200')) then g[i].r200 = rad[hid[i]]
     if (TAG_EXIST(g, 'm200') and TAG_EXIST(h, 'm200')) then g[i].m200 = h[hid[i]].m200
     if (TAG_EXIST(g, 'r200') and TAG_EXIST(h, 'r200')) then g[i].r200 = h[hid[i]].r200
     if KEYWORD_SET(add_centrals) then begin
        if (g[i].rhalo lt 0.001 and g[i].amag(mag_ind) lt h[hid[i]].lbcg) then begin
           g[i].central = 1
           h[hid[i]].lbcg = g[i].amag(mag_ind)
        endif
     endif
     h[hid[i]].ngals++
     h[hid[i]].lumtot += g[i].amag[mag_ind]
     if g[i].central then h[hid[i]].lbcg = g[i].amag[mag_ind]
     if Mr[i] le -18 then h[hid[i]].n18++
     if Mr[i] le -19 then h[hid[i]].n19++
     if Mr[i] le -20 then begin
       h[hid[i]].n20++
       h[hid[i]].lum20 += g[i].amag[mag_ind]
     endif
     if Mr[i] le -21 then h[hid[i]].n21++
     if Mr[i] le -22 then h[hid[i]].n22++
  endif
endfor

print, "# of elements in h: "
help, h

end
