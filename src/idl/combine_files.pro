pro combine_files, path, base, nfiles

gfile = path+'01/idl/'+base+'.01_galaxies.fit'
hfile = path+'01/idl/'+base+'.01_halos.fit'
dfile = path+'01/idl/'+base+'.01_des.fit'
deepfile = path+'01/idl/'+base+'.01_deep.fit'
vfile = path+'01/idl/'+base+'.01_vista.fit'
jfile = path+'01/idl/'+base+'.01_johnson.fit'
afile = path+'01/idl/'+base+'.01_alhambra.fit'

g = mrdfits(gfile, 1)
h = mrdfits(hfile, 1)
d = mrdfits(dfile, 1)
deep = mrdfits(deepfile, 1)
v = mrdfits(vfile, 1)
j = mrdfits(jfile, 1)
a = mrdfits(afile, 1)

for i = 2L, nfiles do begin
  num = '0'+strcompress(string(i),/remove_all)
  num = strmid(num, strlen(num)-2,2)
  
  gfile = path+num+'/idl/'+base+'.'+num+'_galaxies.fit'
  hfile = path+num+'/idl/'+base+'.'+num+'_halos.fit'
  dfile = path+num+'/idl/'+base+'.'+num+'_des.fit'
  deepfile = path+num+'/idl/'+base+'.'+num+'_deep.fit'
  vfile = path+num+'/idl/'+base+'.'+num+'_vista.fit'
  jfile = path+num+'/idl/'+base+'.'+num+'_johnson.fit'
  afile = path+num+'/idl/'+base+'.'+num+'_alhambra.fit'

  tg = mrdfits(gfile, 1)
  th = mrdfits(hfile, 1)
  td = mrdfits(dfile, 1)
  tdeep = mrdfits(deepfile, 1)
  tv = mrdfits(vfile, 1)
  tj = mrdfits(jfile, 1)
  av = mrdfits(afile, 1)

  g = [g, tg]
  h = [h,th]
  d = [d,td]
  deep = [deep,tdeep]
  v = [v,tv]
  j = [j,tj]
  a = [a,ta]

endfor

;;;now we correct the ngals and HOD information
h(*).n19 = 0
h(*).n20 = 0
h(*).n21 = 0
h(*).ngals = 0
h(*).lumtot = 0.
h(*).lum20 = 0.
h(*).lcent = 0.

;;;make a halo lookup table
print, 'Making Lookup Table...'
nht = max(h.haloid)
nh = N_ELEMENTS(h)
ng = N_ELEMENTS(g)
h_table = lonarr(nht+1)
h(*).edge = 0
h(*).edge = 0
for i = 0L, nh - 1 do begin
  h_table(h(i).haloid) = i
endfor

;;;transfer data to halo table
print, 'Transfering Data to Halo Files....'
for i = 0L, ng - 1 do begin
  hind = h_table(g(i).haloid)
  if (g(i).rhalo gt g(i).r200) then continue
  h(hind).ngals++
  h(hind).lumtot += g(i).amag(2)
  if (g(i).amag(2) lt -19.) then h(hind).n19++
  if (g(i).amag(2) lt -20.) then h(hind).n20++
  if (g(i).amag(2) lt -21.) then h(hind).n21++
  if (g(i).central) then h(hind).lcent = g(i).amag(2)
  if (g(i).amag(2) lt -20.) then h(hind).lum20 += g(i).amag(2)
endfor

;;;transfer back to galaxy table
print, 'Transfering back to Galaxies....'
for i = 0L, ng - 1 do begin
  hind = h_table(g(i).haloid)
  if (g(i).rhalo gt g(i).r200) then continue
  g(i).ngals = h(hind).ngals
endfor

;;;write out our files
print, 'Writing out files....'
spawn, 'mkdir '+path+'/idl/'
gfile = path+'/idl/'+base+'_galaxies.fit'
hfile = path+'/idl/'+base+'_halos.fit'
dfile = path+'/idl/'+base+'_des.fit'
deepfile = path+'/idl/'+base+'_deep.fit'
vfile = path+'/idl/'+base+'_vista.fit'
jfile = path+'/idl/'+base+'_johnson.fit'
afile = path+'/idl/'+base+'_alhambra.fit'

mwrfits, g, gfile, /create
mwrfits, h, hfile, /create
mwrfits, d, dfile, /create
mwrfits, deep, deepfile, /create
mwrfits, v, vfile, /create
mwrfits, j, jfile, /create
mwrfits, a, afile, /create

end
