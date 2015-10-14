;pro make_galaxy_ascii_table

for i = 0L, 6 do begin
  file = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/cluster_comp_v1.07_fixedID'+strcompress(string(i),/remove_all)+'.fits'
  g_tmp = mrdfits(file,1)
  if (i eq 0) then g = g_tmp
  if (i gt 0) then g = [g, g_tmp]
  g_tmp = 1
endfor

;ii = where(g.photoz gt 0.1 and g.photoz le 0.3)
;g = g[ii]

h = mrdfits('~/ki01/projects/addgals/PO_oct/v1.07_z/idl/hv_POv1.07_z_oct_halos.fit',1)
ng = N_ELEMENTS(g)
nh = N_ELEMENTS(h)
incl = bytarr(nh)

;;;make lookup table
print, 'Making first lookup table'
table = lonarr(max(h.haloid)+1)
for i = 0L, nh - 1 do begin
  table(h(i).haloid) = i
endfor

;;;find halos in the sample
print, 'Getting Halo Sample'
incl(*) = 0
for i = 0L, ng - 1 do begin
  if (g(i).haloid ge 0) then hid = table(g(i).haloid)
  incl(hid) = 1
endfor

ii = where(incl eq 1)
h = h[ii]
srt = reverse(sort(h.m200))
h = h[srt]

;;;remake the lookup table
table = lonarr(max(h.haloid)+1)
nh = N_ELEMENTS(h)
for i = 0L, nh - 1 do begin
  table(h(i).haloid) = i
endfor

g_str = create_struct('rank', 0L, 'GalaxyID', 0L)
g_ascii = replicate(g_str, ng)

print, 'getting ranks'
for i = 0L, ng - 1 do begin
  rank = -1
  if (g(i).haloid ge 0) then rank = table(g(i).haloid)+1
  g_ascii(i).rank = rank
endfor
g_ascii.GalaxyID = g.id
srt = sort(g_ascii.rank)
g_ascii = g_ascii[srt]

;;;write the ascii files
out_halo = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/halo2.ascii'
out_galaxy = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/galaxy2.ascii'

print, 'writting files'
openw,1,out_halo
for i = 0L, nh - 1 do begin
  printf,1,i+1, h[i].ra, h[i].dec, h[i].z, h[i].m200, h[i].ngals
endfor
close,1

openw,1,out_galaxy
for i = 0L, ng - 1 do begin
  printf, 1, g_ascii[i].rank, g_ascii[i].galaxyid
endfor
close,1

;;;look at the actual # of galaxies found in each halo
ngals_found = lonarr(nh)
ngals_found(*) = 0
for i = 0L, ng - 1 do begin
  if (g(i).haloid lt 0) then continue
  hid = table(g(i).haloid)
  ngals_found(hid)++
endfor  

end
