g = mrdfits('~/ki01/projects/addgals/PO_oct/v1.07_z/idl/g_tmp.fit',1)
g_dim = mrdfits('~/ki01/projects/addgals/PO_oct/v1.07_z/idl/g_dim_tmp.fit',1)

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

;;;kludge to fiz galaxy id
for i = 0L, ng - 1 do begin
  g_tot.id = i
endfor

print, "Splitting and writing out the files."
z_limit = [0.1, 0.3, 0.4, 0.6, 0.7, 0.9, 1.1, 1.3]
nbins = N_ELEMENTS(z_limit) - 1

for i = 0L, nbins - 1 do begin
  print, i
  ii = where(g_tot.photoz ge z_limit(i) and g_tot.photoz lt z_limit(i+1))
  outfile = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/cluster_comp_v1.07_fixedID_'+strcompress(string(i),/remove_all)+'.fits'

  g_tmp = g_tot[ii]
  REMOVE_TAGS, g_tmp, ['tmag', 'omag', 'omagerr'], g_new0
  add_tag, g_new0, 'tmag', fltarr(5), g_new1
  g_new0 = 0
  add_tag, g_new1, 'omag', fltarr(5), g_new2
  g_new1 = 0
  add_tag, g_new2, 'omagerr', fltarr(5), g_new
  g_new2 = 0
  g_new.tmag(0:3) = g_tmp.tmag(0:3)
  g_new.tmag(4) = g_tmp.tmag(5)
  g_new.omag(0:3) = g_tmp.omag(0:3)
  g_new.omag(4) = g_tmp.omag(5)
  g_new.omagerr(0:3) = g_tmp.omagerr(0:3)
  g_new.omagerr(4) = g_tmp.omagerr(5)

  ;;;add amags
  filters = ['DES_g.par', 'DES_r.par', 'DES_i.par', 'DES_z.par', 'DES_Y.par']
  kcorrect, g_new.omag, g_new.omagerr, g_new.photoz, kcorrect, $
            band_shift = 0.1,$
            coeffs=coeffs,/magnitude,/stddev,filterlist=filters
  dm = lf_distmod(g_new.photoz)
  add_tag, g_new, 'amag', fltarr(5), gn
  for k = 0, 4 do gn.amag(k) = gn(*).omag(k) - dm(*) - kcorrect(k,*)
  g_new = gn
  gn = 0

  mwrfits, g_new, outfile, /create
  g_new = 0
endfor

end

