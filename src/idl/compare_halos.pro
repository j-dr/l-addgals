bdir = '/nfs/slac/g/ki/ki01/mbusha/projects/addgals/Carmen/'

zmin = 0.1
zmax = 0.33
mmin = 5e13

dir1 = bdir + 'LC/v1.08a/v16/catalog/'
dir2 = bdir + 'LC/v1.08a/v16/DES_catalog/'
dir3 = bdir + 'LC/v1.08a/v21/DES_catalog/'

file1 = dir1 + 'halos_24.fit'
file2 = dir2 + 'DES_Mock_halos_24.fit'
file3 = dir3 + 'DES_Mock_halos_24.fit'

h1 = mrdfits(file1,1)
h2 = mrdfits(file2,1)
h3 = mrdfits(file3,1)

;;;cut to limits
ii1 = where(h1.m200 gt mmin and h1.z gt zmin and h1.z lt zmax, count1)
ii2 = where(h2.m200 gt mmin and h2.z gt zmin and h2.z lt zmax, count2)
ii3 = where(h3.m200 gt mmin and h3.z gt zmin and h3.z lt zmax, count3)

print, "Number of centrals in our samples:  ", count1, count2, count3

;;;histograms of mass
window, 0
bsize = 0.03
plothist, alog10(h1[ii1].m200), bin = bsize, /ylog, xtitle = 'log(M200)', $
          yrange = [0.1,1e2]
plothist, alog10(h2[ii2].m200), bin = bsize, /over, color = !red, linestyle=2
plothist, alog10(h3[ii3].m200), bin = bsize, /over, color = !blue, linestyle=1

window, 1
zbsize = 0.01
plothist, h1[ii1].z, bin = zbsize, xtitle = 'z'
plothist, h2[ii2].z, bin = zbsize, /over, color = !red, linestyle = 2
plothist, h3[ii3].z, bin = zbsize, /over, color = !blue, linestyle = 1

end
