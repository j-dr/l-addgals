pro compare_wp, dir, sfx=sfx, linestyle=linestyle, label=label

; options are typicall all, red, and blue
if NOT KEYWORD_SET(sfx) then sfx = 'all'

if NOT KEYWORD_SET(linestyle) then linestyle = lindgen(N_ELEMENTS(dir))

; set what we're going to plot
nsfx = N_ELEMENTS(sfx)
mag = ['19', '20', '21']
;mag = ['19', '20', '21']

for i = 0L, N_ELEMENTS(mag) - 1 do begin

  ; setup our i/o files
  tail='/wp_mt'+mag[i]+'_'+sfx+'.dat'

  ; setup a blank plot
  plot, [1,1], [1,1], /xlog, /ylog, xrange = [0.1,40], yrange = [1,2e3], $
	xtitle = 'r!Dp!N [Mpc/h]', ytitle = 'w!Dp!N(r!Dp!N)', $
	title = 'M!Dr!N < -'+mag[i], /nodat, /xsty, /ysty

  ; loop through the sfx's to add
  for idir = 0L, N_ELEMENTS(dir) - 1 do begin
    f = dir[idir]+tail
      for j = 0L, nsfx - 1 do begin

      ; set the plot color (black, red, or blue)
      pcolor = !P.COLOR
      if (sfx[j] eq 'red') then pcolor = !red
      if (sfx[j] eq 'blue') then pcolor = !blue

      ; read the data
      data = read_wp(f[j])

      ; add to the plot
      oploterror, data.r, data.xi, data.xi_err, color=pcolor, errcolor=pcolor, $
	linestyle = linestyle[idir]
    endfor
  endfor
  if KEYWORD_SET(label) then legend, label, linestyle=linestyle, /right
endfor


end
    
  
