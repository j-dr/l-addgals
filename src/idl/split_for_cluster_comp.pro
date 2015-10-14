pro split_for_cluster_comp, file

bright = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/hv_POv1.07_z_des_5year.fit'
halos = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/hv_POv1.07_z_oct_halos.fit'
g = mrdfits(bright,1)
ii = where(g.ra gt 0. and g.ra lt 10.)
g = g[ii]
h = mrdfits(halos,1)
dim1 = '~/ki01/data/des/official_catalogs/v1.06/hv_POv1.06b_oct.01_observed_5year_dim.fit'
dim2 = '~/ki01/data/des/official_catalogs/v1.06/hv_POv1.06b_oct.02_observed_5year_dim.fit'
g_dim1 = mrdfits(dim1, 1)
g_dim2 = mrdfits(dim2, 1)


;;;rescale r200 data for halos
table = lindgen(max(h.haloid)+1)
for i = 0L, N_ELEMENTS(h) - 1 do begin
  table(h(i).haloid) = i
  a = 1./(1+h[i].z)
  h[i].r200 *= a
endfor
for i = 0L, N_ELEMENTS(g)-1 do begin
  hind = table[g[i].haloid]
  g[i].r200 = h[hind].r200
endfor

table = 0

;g_dim1.ra += 1.
;g_dim2.ra += 1.
;ii = where(g_dim1.ra lt 1.)
;g_dim1[ii].ra += 10.

ii = where(g_dim1.omag(2) lt 25. and g_dim1.photoz lt 1.3)
g_dim1 = g_dim1[ii]
ii = where(g_dim2.omag(2) lt 25. and g_dim2.photoz lt 1.3)
g_dim2 = g_dim2[ii]

g_dim = [g_dim1, g_dim2]
g_dim1 = 0
g_dim2 = 0
g_dim.ra *= (1+randomn(seed,N_ELEMENTS(g_dim))*0.01)
ii = where(g_dim.ra gt 0. and g_dim.ra lt 10.)
g_dim = g_dim[ii]

;;;find the nearest halos for the dim galaxies
print, "Finding nearest halos for dim galaxies"
bins = 64
MakeLinkedList,x=h.halopx, y=h.halopy, z=h.halopz, llbins=bins,llhoc=llhoc,$
               llist=llist,boxsize=3000.
ndim = N_ELEMENTS(g_dim)
;;;convert to cartesian relatively effeciently -- redshift 
z_table = generate_z_of_r_table(0.3, 0.7, zmax = 1.5)
r = fltarr(ndim)
;for i = 0L, ndim - 1 do begin
;  r(i) = r_of_z(g_dim(i).z, z_table)
;endfor
r = r_of_z(g_dim.z, z_table)
z = sin(g_dim.dec*!PI/180.)*r
xy = cos(g_dim.dec*!PI/180.)*r
y = sin(g_dim.ra*!PI/180.)*xy
x = cos(g_dim.ra*!PI/180.)*xy
scalefac = bins/3000.
for i = 0L, ndim - 1 do begin
  ix = floor(x(i)*scalefac)
  iy = floor(y(i)*scalefac)
  iz = floor(z(i)*scalefac)
  r2_closest = 1e10
  ih_closest = -1
  r2 = 1e10
  for iix = ix-1,ix+1 do begin
    if (iix lt 0 or iix ge bins) then continue
    for iiy = iy-1,iy+1 do begin
      if (iiy lt 0 or iiy ge bins) then continue
      for iiz = iz-1,iz+1 do begin
        if (iiz lt 0 or iiz ge bins) then continue
        ih = LLHoc(iix,iiy,iiz)
        while (ih ge 0) do begin
          dx = x(i) - h[ih].halopx
          dy = y(i) - h[ih].halopy
          dz = z(i) - h[ih].halopz
          r2 = dx*dx + dy*dy + dz*dz
          if (r2 lt r2_closest) then begin
            r2_closest = r2
            ih_closest = ih
          endif
          ih = llist(ih)
        endwhile
      endfor
    endfor
  endfor
  g_dim[i].rhalo = sqrt(r2_closest)
  g_dim[i].haloid = h[ih_closest].haloid
  g_dim[i].r200 = h[ih_closest].r200
endfor

mwrfits, g, '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/g_tmp.fit', /create
mwrfits, g_dim, '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/g_dim_tmp.fit', /create

;;;reset the hid stuff
print, "Correcting halo occupancy"
g_tot = [g, g_dim]
g = 0
g_dim = 0
ng = N_ELEMENTS(g_tot)
ii = where(g_tot.rhalo gt g_tot.r200)
g_tot[ii].haloid = -1.
g_tot[ii].r200 = -1.
g_tot[ii].rhalo = -1.

;;;remove the Z band
print, "Removing z-band"
REMOVE_TAGS, g_tot, ['tmag', 'omag', 'omagerr'], g_new
add_tag, g_new, 'tmag', fltarr(5), g_new1
add_tag, g_new1, 'omag', fltarr(5), g_new2
g_new1 = 0
add_tag, g_new2, 'omagerr', fltarr(5), g_new3
g_new2 = 0
g_new = g_new3
g_new3 = 0
g_new.tmag(0:3) = g.tmag(0:3)
g_new.tmag(4) = g.tmag(5)
g_new.omag(0:3) = g.omag(0:3)
g_new.omag(4) = g.omag(5)
g_new.omagerr(0:3) = g.omagerr(0:3)
g_new.omagerr(4) = g.omagerr(5)

print, "Splitting and writing out the files."
z_limit = [0.1, 0.3, 0.4, 0.6, 0.7, 0.9, 1.1, 1.3]
nbins = N_ELEMENTS(z_limit) - 1
for i = 0L, nbins - 1 do begin
  ii = where(g_new.photoz ge z_limit(i) and g_new.photoz lt z_limit(i+1))
  outfile = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/cluster_comp_v1.07_'+strcompress(string(i),/remove_all)+'.fits'
  mwrfits, g_new[ii], outfile, /create
endfor

end

