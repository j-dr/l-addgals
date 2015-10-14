pro compare_centrals_raw_files

bdir = '/nfs/slac/g/ki/ki01/mbusha/projects/addgals/Carmen/'

zmin = 0.1
zmax = 0.33
mmin = 5e13
magmin = -19.

dir1 = bdir + 'LC/v1.08a/v16/'
dir2 = bdir + 'LC/v1.08a/v22/'

gfile = '/hv_output/gal_ginfo1.dat'
hfile = '/hv_output/gal_hinfo.dat'

nfiles1 = 3
nfiles2 = 1

read_files, dir1, gfile, hfile, nfiles1, sed1, mr1, ra1, dec1, z1, cent1, m2001
read_files, dir2, gfile, hfile, nfiles2, sed1, mr2, ra2, dec2, z2, cent2, m2002

;;;cut to limits
cc1 = where(cent1 and z1 gt zmin and z1 lt zmax and mr1 lt magmin, count1)
cc2 = where(cent2 and z2 gt zmin and z2 lt zmax and mr2 lt magmin, count2)
print, "Number of centrals in our samples:  ", count1, count2

gg1 = where(z1 gt zmin and z1 lt zmax and mr1 lt magmin, gcount1)
gg2 = where(z2 gt zmin and z1 lt zmax and mr2 lt magmin, gcount2)
print, "Number of galaxies in our samples:  ", gcount1, gcount2

;;;histograms of mass
window, 0
bsize = 0.03
plothist, alog10(m2001[cc1]), bin = bsize, /ylog, xtitle = 'log(M200)', $
          yrange = [0.1,1e2], thick = 3
plothist, alog10(m2002[cc2]), bin = bsize, /over, color = !red, linestyle=2, $
          thick = 3

window, 1
mbinsize = 0.1
plothist, mr1[gg1], bin = mbinsize, xtitle = 'M!Dr', /ylog, yrange=[0.1,1e4]
plothist, mr1[cc1], bin = mbinsize, xtitle = 'M!Dr', /over, color = !red
plothist, mr2[gg1], bin = mbinsize, xtitle='M!Dr',/over,linestyle=2, thick=3
plothist, mr2[cc2], bin = mbinsize, xtitle='M!Dr',/over,color=!red,linestyle=2

window,2
plot, m2001[cc1], mr1[cc1], psym = 1, /xlog, xtitle = 'M!D200', ytitle='M!Dr'
oplot, m2002[cc2], mr2[cc2], psym = 1, color = !red

mstar1 = -20.73
mstar2 = -20.41
lum1 = 10.^(-0.4*(mr1 - mstar1))
lum2 = 10.^(-0.4*(mr2 - mstar2))
plot, m2001[cc1], lum1[cc1], psym = 1, /xlog, /ylog, xtitle = 'M!D200', ytitle='L',$
      xrange = [5e13, 1e15], /xsty, yrange = [1,20], /ysty
oplot, m2002[cc2], lum2[cc2]*1.2, psym = 1, color = !red

end

pro read_files, bdir, gfile, hfile, nfiles, sed, mr, ra, dec, z, cent, m200

for i = 0, nfiles - 1 do begin
   num = strcompress(string(i), /remove_all)
   rdfloat, bdir+num+gfile, tsed, tmr, tra, tdec, tz, tcent
   rdfloat, bdir+num+hfile, id, tm200
   if (i eq 0) then begin
      mr = tmr
      ra = tra
      dec = tdec
      z = tz
      cent = tcent
      m200 = tm200
   endif else begin
      mr = [mr, tmr]
      ra = [ra, tra]
      dec = [dec, tdec]
      z = [z, tz]
      cent = [cent, tcent]
      m200 = [m200, tm200]
   endelse
endfor

cent = floor(cent)

end
