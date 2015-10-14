pro cheap_lf, g, g_dim

z1 = 0.09
z2 = 0.11

ind = where(g.z gt z1 and g.z lt z2)
ind_dim = where(g_dim.z gt z1 and g_dim.z lt z2)

r1= lookbackdist(z1, 0.3, 0.7)
r2= lookbackdist(z2, 0.3, 0.7)
vol = 4./3*!PI*(r2^3-r1^3)/8

fit_lf, [g(ind).amag(1), g_dim(ind_dim).amag(1)], box=vol^(1./3), gmagbin1, gdensity1
fit_lf, [g(ind).amag(2), g_dim(ind_dim).amag(2)], box=vol^(1./3), gmagbin2, gdensity2
fit_lf, [g(ind).amag(3), g_dim(ind_dim).amag(3)], box=vol^(1./3), gmagbin3, gdensity3

xr = [-23, -18]
yr = [1e-6, 1e-1]

plot, gmagbin1, gdensity1, /ylog, yrange=yr, xrange=xr, /ysty, xtitle='M-5logh', ytitle=textoidl('\phi [h!u3!nMpc!u-3!n]'), /xsty, psym = 2, /nodat
oplot, gmagbin1, gdensity1, psym = 2, color = !green
oplot, gmagbin2, gdensity2, psym = 2, color = !red
oplot, gmagbin3, gdensity3, psym = 2, color = !blue

blan = [2.18e-2,-0.89, -19.39]
oplot, gmagbin1, schechter_mag(gmagbin1, blan), color=!green
blan = [1.49e-2,-1.05, -20.44]
oplot, gmagbin2, schechter_mag(gmagbin2, blan), Color=!red
blan = [1.47e-2,-1.00, -20.82]
oplot, gmagbin3, schechter_mag(gmagbin3, blan), color=!blue

return
end
