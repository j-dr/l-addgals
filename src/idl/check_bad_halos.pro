pro check_bad_halos, h1, h2

mlim = 6.8e13
slope = alog10(400)/alog10(1e15/2e12)
yint = -alog10(2e12)*slope

ii1 = where(h1.z gt 0.6 and h1.z lt 0.7)
ii2 = where(h2.z gt 0.6 and h2.z lt 0.7)

plot, h1[ii1].m200, h1[ii1].n18, psym = 1, /xlog, /ylog, yrange = [1,1e4], xtitle = 'M!D200', ytitle = 'n!D18'
asdf
oplot, h2[ii2].m200, h2[ii2].n18, psym = 1, color = !red

oplot, [2e12,1e15], [1,400], color = !blue
xpts = [2e12,1e15]
ypts = 10.^(slope*alog10(xpts) + yint)
oplot, xpts, ypts, color = !blue
oplot, [1,1]*mlim, [1,1e4], color = !blue

bad1 = ii1[where(h1[ii1].m200 gt mlim and alog10(h1[ii1].n18) gt yint+slope*alog10(h1[ii1].m200), count1)]
bad2 = ii2[where(h2[ii2].m200 gt mlim and alog10(h2[ii2].n18) gt yint+slope*alog10(h2[ii2].m200), count2)]

print, count1, count2

match, h1[bad1].haloid, h2.haloid, a1, a2
match, h1.haloid, h2[bad2].haloid, b1, b2

print, mean(float(h1[b1].n18)/h2[bad2[b2]].n18)
range, float(h1[b1].n18)/h2[bad2[b2]].n18


end
