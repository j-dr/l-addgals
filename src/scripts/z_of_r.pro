function z_of_r, r, table

npts = N_ELEMENTS(table)
nz = N_ELEMENTS(r)
zred = fltarr(nz)
print, 'npts: ', npts
if (nz eq 1) then begin
   for i = 1L, npts - 1 do begin
      if (table(i).r gt r) then break
   endfor
   print, i
   slope = (table(i).z - table(i-1).z)/(table(i).r-table(i-1).r)
   zred = table(i-1).z + slope*(r-table(i-1).r)
endif else begin
   for i = 1L, npts - 1 do begin
      ii = where(r ge table(i-1).r and r lt table(i).r, count)
      if (count eq 0) then continue
      slope = (table(i).z - table(i-1).z)/(table(i).r-table(i-1).r)
      zred(ii) = table(i-1).z + slope*(r(ii)-table(i-1).r)
   endfor
endelse

return, zred


end
