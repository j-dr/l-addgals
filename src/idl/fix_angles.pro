pro fix_angles, g

dec  = g.ra
ra = g.dec
ii = where(ra lt 0)
if (ii[0] ne -1) then ra[ii] = ra[ii]+360.
g.ra = ra
g.dec = dec


dec  = g.halora
ra = g.halodec
if (ii[0] ne -1) then ra[ii] = ra[ii]+360.
g.halora = ra
g.halodec = dec
end

