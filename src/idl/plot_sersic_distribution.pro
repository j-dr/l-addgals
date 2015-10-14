function plot_sersic_distribution, sersic

carr = [!yellow, !green, !red, !blue, !purple]

bin = 0.1
pplothist, sersic.sersic_n[4], bin=bin, xrange = [-1, 7], xtitle = 'sersic index'
for i = 0, 4 do begin
  pplothist, sersic.sersic_n[i], bin=bin, /over, color = carr[i]
  ka = where(sersic.sersic_n[i] eq max(sersic.sersic_n), count)
  print, "Fraction with maximum sersic index: ", float(count)/float(N_ELEMENTS(sersic))
endfor
legend, ['u', 'g', 'r', 'i', 'z'], linestyle = 0, color=carr

bin = 0.3
pplothist, sersic.sersic_r0[0], bin=bin, xtitle = 'r0', /ylog, yrange = [1,1e5]
for i = 0, 4 do begin
  pplothist, sersic.sersic_r0[i], bin=bin, /over, color = carr[i]
  ka = where(sersic.sersic_r0[i] eq max(sersic.sersic_r0), count)
  print, "Fraction with maximum r0: ", float(count)/float(N_ELEMENTS(sersic))
endfor
legend, ['u', 'g', 'r', 'i', 'z'], linestyle = 0, color=carr, /right

bin = 0.3
pplothist, sersic.sersic_r50[0], bin=bin, xtitle = 'r50', /ylog, yrange = [1,1e5]
for i = 0, 4 do begin
  pplothist, sersic.sersic_r50[i], bin=bin, /over, color = carr[i]
  ka = where(sersic.sersic_r50[i] eq max(sersic.sersic_r50), count)
  print, "Fraction with maximum r50: ", float(count)/float(N_ELEMENTS(sersic))
endfor
legend, ['u', 'g', 'r', 'i', 'z'], linestyle = 0, color=carr, /right

bin = 5
pplothist, sersic.sersic_r90[0], bin=bin, xtitle = 'r90', /ylog, yrange = [1,1e5]
for i = 0, 4 do begin
  pplothist, sersic.sersic_r90[i], bin=bin, /over, color = carr[i]
  ka = where(sersic.sersic_r90[i] eq max(sersic.sersic_r90), count)
  print, "Fraction with maximum r90: ", float(count)/float(N_ELEMENTS(sersic))
endfor
legend, ['u', 'g', 'r', 'i', 'z'], linestyle = 0, color=carr, /right

; select objects that have at least one decent measurement in gri
max_sersic = max(sersic.sersic_n)
nbad = 0L
good_measure = intarr(N_ELEMENTS(sersic))
for i = 0L, N_ELEMENTS(sersic) - 1 do begin
  ka = where(sersic[i].sersic_n[1:4] lt max_sersic, count)
  if (count eq 0) then begin
    nbad++
    good_measure[i] = 1
  endif
endfor
print, "Fraction of galaxies with no good sersic measurement in griz: ", float(nbad)/float(N_ELEMENTS(sersic))
ka = where(good_measure eq 1)

return, ka

end
