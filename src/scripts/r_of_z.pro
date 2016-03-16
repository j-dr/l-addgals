function r_of_z, z, table

npts = N_ELEMENTS(table)

;for i = 1L, npts - 1 do begin
;  if (table(i).z gt z) then break
;endfor
;slope = (table(i).r - table(i-1).r)/(table(i).z-table(i-1).z)
;rofz = table(i-1).r + slope*(z-table(i-1).z)

rofz = fltarr(N_ELEMENTS(z))
rofz(*) = -1
for i = 0L, npts - 2 do begin
  ii = where(z ge table(i).z and z lt table(i+1).z)
  if (ii(0) eq -1) then continue
  slope = (table(i+1).r - table(i).r)/(table(i+1).z - table(i).z)
  rofz(ii) = table(i).r + slope*(z(ii)-table(i).z)
endfor

if (N_ELEMENTS(z) eq 1) then rofz = rofz[0]

return, rofz

end
