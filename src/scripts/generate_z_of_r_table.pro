function generate_z_of_r_table, omegam, omegal, zmax=zmax, npts=npts

if not KEYWORD_SET(npts) then npts = 1000
;if not (arg_present(zmax)) then zmax = 2.0
if not (KEYWORD_SET(zmax)) then zmax = 2.0
c = 2.9979e5
da = (1.0 - (1.0/(zmax+1)))/npts
table = create_struct('r', 0.0, 'z', 0.0)
z_of_r_table = replicate(table, npts)
Thisa = 1.0
z_of_r_table(0).z = 1.0/Thisa - 1.
z_of_r_table(0).r = 0.0
for i = 1L, npts - 1 do begin
  Thisa = 1. - da*float(i)
  ThisH = 100.*sqrt(omegam/Thisa^3 + omegal)
  z_of_r_table(i).z = 1./Thisa - 1
  z_of_r_table(i).r = z_of_r_table(i-1).r + 1./(ThisH*Thisa*Thisa)*da*c
endfor

return, z_of_r_table

end
