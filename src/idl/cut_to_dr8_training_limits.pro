function cut_to_dr8_training_limits, g, ii=ii

ii1 = where(g.omag[2] lt 21.8, count1)

for i = 22.1, 32.1, 0.1 do begin
  ii2 = where(g[ii1].omag[0] lt i and g[ii1].omag[1] lt 29 and $
              g[ii1].omag[3] lt 29 and g[ii1].omag[4] lt 29 and $
	      g[ii1].omag[0] gt 0 and g[ii1].omag[1] gt 0 and $
              g[ii1].omag[3] gt 0 and g[ii1].omag[4] gt 0, count2)
  if (float(count2)/float(count1) gt 0.8) then break
endfor
print, 'Using u-band cut ', i


ii = ii1[ii2]
gp = g[ii]

return, gp

end
