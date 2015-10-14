
;yint = 41.8
;slope = 4.0
yint = -25.0-0.2
slope = -5.2

mstar = -20.44

z = 0.01*(findgen(50)+1)
mmin_of_z = yint+alog10(z)*slope
mstar_of_z = mstar - 1.3*(z-0.1)
;mstar_of_z(*) = mstar
;dm = lf_distmod(z)
;amagmin_of_z = mmin_of_z - dm
;lmin_of_z = 10.^(-0.4*(amagmin_of_z-mstar_of_z))
lmin_of_z = 10.^(-0.4*(mmin_of_z-mstar_of_z))
bbstar = 0.85 + 0.15*lmin_of_z - 0.04*(mmin_of_z-mstar_of_z)

plot, z, bbstar, xtitle='z', ytitle='b/b!D*'
;plot,mmin_of_z,bbstar,xtitle='M!Dr',ytitle='b/b!D*',xrange=[-17,-23],/xsty,yrange=[0,2.4],/ysty


end
