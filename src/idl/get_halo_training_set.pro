;pro get_halo_training_set, path, catname

dir = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/Buzzard/Catalog_v1.0/'
catname = 'Buzzard_v1.0'
skip_read = 1

num = ['86', '114', '115', '147']
gfiles = dir+'/truth_rotated/'+catname+'_truth_des.'+num+'.fit'
goutfile = dir+'/truth_rotated/'+catname+'_truth_des_redmapper_sepc.fit'
hfile = dir+'/halos/'+catname+'_halos_rotated.4.fit'
houtfile = dir+'/halos/'+catname+'_halos_redmapper_sepc.fit'

if (skip_read eq 0) then begin
g = mrdfits(gfiles[0], 1)
for i = 1, 3 do begin
  tg = mrdfits(gfiles[i],1)
  g = [g, tg]
endfor

h = mrdfits(hfile,1)
cc = where(h.haloid eq h.host_haloid and h.m200 gt 5e13)
hc1 = h[cc]
endif

;;;cut at the galaxy magnitude limit
gc1 = g[where(g.omag[1] lt 20 and g.central)]

;;;match our galaxies to our halos
match, gc1.haloid, hc1.haloid, a, b
gc = gc1[a]
hc = hc1[b]
nh = N_ELEMENTS(hc)
keep = intarr(nh)
hc.ra = gc.ra
hc.dec = gc.dec

;;;select galaxies
hist = histogram(hc.z, bin = 0.01, min = 0.0, max = 0.7, reverse_indices=ri)
nkeep_bin = 30
for i = 0, N_ELEMENTS(hist)-1 do begin
  if (hist[i] le 0) then continue
  ind = ri[ri[i]+lindgen(hist[i])]
  if (hist[i] le nkeep_bin) then ind2 = ind else $
    ind2 = get_unique_sample_number(ind, nkeep_bin)
  keep[ind2] = 1 
endfor
ka = where(keep)
hout1 = hc[ka]
gout1 = gc[ka]

;;;cut at the galaxy magnitude limit
theta = (90.-g.dec)*!PI/180.
phi = g.ra*!PI/180.
ang2pix_ring, 16, theta, phi, ip
gc1 = g[where(g.omag[1] lt 22 and g.central and (ip eq 2012 or ip eq 2077))]

;;;match our galaxies to our halos
match, gc1.haloid, hc1.haloid, a, b
gc = gc1[a]
hc = hc1[b]
nh = N_ELEMENTS(hc)
keep = intarr(nh)
hc.ra = gc.ra
hc.dec = gc.dec

hist = histogram(hc.z, bin = 0.01, min = 0.7, max = 1.1, reverse_indices=ri)
nkeep_bin = 10
for i = 0, N_ELEMENTS(hist)-1 do begin
  if (hist[i] le 0) then continue
  ind = ri[ri[i]+lindgen(hist[i])]
  if (hist[i] le nkeep_bin) then ind2 = ind else $
    ind2 = get_unique_sample_number(ind, nkeep_bin)
  keep[ind2] = 1
endfor
ka = where(keep)
hout2 = hc[ka]
gout2 = gc[ka]

hout = [hout1, hout2]
gout = [gout1, gout2]

mwrfits, hout, houtfile, /create
mwrfits, gout, goutfile, /create

end
