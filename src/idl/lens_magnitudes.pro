pro lens_magnitudes, g, shear, nmag

;;;add the lensing
if not (tag_exist(g, 'lmag')) then begin
  add_tags, g, ['lmag'], ['fltarr('+string(nmag)+')'], g2
  g = g2
  g2 = 0
endif
g.lmag = 0.
mu = 1./(shear.a11*shear.a00 - shear.a01*shear.a10)
for ig = 0L, N_ELEMENTS(shear) - 1 do begin
  tind = shear[ig].index
  if (g[tind].lmag[0] ne 0.) then begin
    g = [g, g[tind]]
    tind = N_ELEMENTS(g) - 1
  endif
  for im = 0, nmag - 1 do $
    g[tind].lmag[im] = g[tind].tmag[im] - 2.5*alog10(mu[ig])
endfor

end

