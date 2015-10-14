bdir = '/nfs/slac/g/ki/ki01/mbusha/projects/addgals/Carmen/'

zmin = 0.1
zmax = 0.33
mmin = 5e13

dir1 = bdir + 'LC/v1.08a/v16/catalog/'
dir2 = bdir + 'LC/v1.08a/v16/DES_catalog/'
dir3 = bdir + 'LC/v1.08a/v21/DES_catalog/'

file1 = dir1 + 'galaxies_21.fit'
;file2 = dir2 + 'DES_Mock_galaxies_21.fit'
file3 = dir3 + 'DES_Mock_galaxies_21.fit'
hfile1 = dir1 + 'halos_24.fit'
hfile2 = dir2 + 'DES_Mock_halos_24.fit'
hfile3 = dir3 + 'DES_Mock_halos_24.fit'

g1 = mrdfits(file1,1)
;g2 = mrdfits(file2,1)
g3 = mrdfits(file3,1)
h1 = mrdfits(hfile1,1)
h2 = mrdfits(hfile2,1)
h3 = mrdfits(hfile3,1)

;;;cut to limits
cc1 = where(g1.central and g1.z gt zmin and g1.z lt zmax, count1)
;cc2 = where(g2.central and g2.z gt zmin and g2.z lt zmax, count2)
cc3 = where(g3.central and g3.z gt zmin and g3.z lt zmax, count3)
print, "Number of centrals in our samples:  ", count1, count3

gg1 = where(g1.z gt zmin and g1.z lt zmax, gcount1)
;gg2 = where(g2.z gt zmin and g2.z lt zmax, gcount2)
gg3 = where(g3.z gt zmin and g3.z lt zmax, gcount3)
print, "Number of galaxies in our samples:  ", gcount1, gcount3

ii1 = where(h1.m200 gt mmin and h1.z gt zmin and h1.z lt zmax, hcount1)
ii2 = where(h2.m200 gt mmin and h2.z gt zmin and h2.z lt zmax, hcount2)
ii3 = where(h3.m200 gt mmin and h3.z gt zmin and h3.z lt zmax, hcount3)

print, "Number of halos in our samples:  ", hcount1, hcount2, hcount3

;;;histograms of mass
window, 0
bsize = 0.03
plothist, alog10(g1[cc1].m200), bin = bsize, /ylog, xtitle = 'log(M200)', $
          yrange = [0.1,1e2], thick = 3
;plothist, alog10(g2[cc2].m200), bin = bsize, /over, color = !red, linestyle=2, $
;          thick = 3
plothist, alog10(g3[cc3].m200), bin = bsize, /over, color = !blue, linestyle=1, $
          thick = 3
plothist, alog10(h1[ii1].m200), bin = bsize, /over, yrange = [0.1,1e2]
;plothist, alog10(h2[ii2].m200), bin = bsize, /over, color = !red, linestyle=2
;plothist, alog10(h3[ii3].m200), bin = bsize, /over, color = !blue, linestyle=1

window, 1
zbsize = 0.01
plothist, g1[cc1].z, bin = zbsize, xtitle = 'z', thick = 3
;plothist, g2[cc2].z, bin = zbsize, /over, color = !red, linestyle = 2, thick=3
plothist, g3[cc3].z, bin = zbsize, /over, color = !blue, linestyle = 1, thick=3
plothist, h1[ii1].z, bin = zbsize, /over
;plothist, h2[ii2].z, bin = zbsize, /over, color = !red, linestyle = 2
;plothist, h3[ii3].z, bin = zbsize, /over, color = !blue, linestyle = 1

mbinsize = 0.1
plothist, g1[gg1].amag[2], bin = mbinsize, xtitle = 'M!Dr', /ylog, yrange=[0.1,1e4]
plothist, g1[cc1].amag[2], bin = mbinsize, xtitle = 'M!Dr', /over, color = !red
plothist, g3[gg1].amag[2], bin = mbinsize, xtitle='M!Dr',/over,linestyle=2, thick=3
plothist, g3[cc3].amag[2]-0.6, bin = mbinsize, xtitle='M!Dr',/over,color=!red,linestyle=2

window,2
plot, g1[cc1].m200, g1[cc1].amag[2], psym = 3, /xlog, xtitle = 'M!D200', ytitle='M!Dr'
oplot, g3[cc3].m200, g3[cc3].amag[2], psym = 3, color = !red

end
