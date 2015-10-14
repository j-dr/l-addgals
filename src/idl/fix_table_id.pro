
first_gid = 1L

for i = 0L, 6 do begin
  file = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/cluster_comp_v1.07_'+strcompress(string(i),/remove_all)+'.fits'
  g = mrdfits(file,1)
  ng = N_ELEMENTS(g)
  g.id = lindgen(ng)+first_gid
  ofile = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/cluster_comp_v1.07_fixedID'+strcompress(string(i),/remove_all)+'.fits'
  mwrfits, g, ofile, /create

  file = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/luster_comp_v1.07_amag_'+strcompress(string(i),/remove_all)+'.fits'
  g = mrdfits(file,1)
  g.id = lindgen(ng)+first_gid
  ofile = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/cluster_comp_v1.07_amag_fixedID'+strcompress(string(i),/remove_all)+'.fits'
  mwrfits, g, ofile, /create

  first_gid += ng
endfor

end
