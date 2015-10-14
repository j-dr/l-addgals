pro make_bcc_halo_catalog, path, inbase, outname, pixnum, zmin=zmin, zmax=zmax, mmin=mmin

if not KEYWORD_SET(zmin) then zmin = 0.0
if not KEYWORD_SET(zmax) then zmax = 10.0
if not KEYWORD_SET(mmin) then mmin = 0.0

nfiles = N_ELEMENTS(pixnum)

h = mrdfits(path+'/'+inbase+'.'+pixnum[0]+'.fit', 1)
range, h.ngals
for i = 1, nfiles - 1 do begin
  print, i, ' ', pixnum[i]
  tfile = path+inbase+'.'+pixnum[i]+'.fit'
;  th = mrdfits(tfile, 1)
  th = read_fits_struct_tags(tfile, ['ngals','lbcg','n18','n19','n20','n21','n22','lumtot','lum20'])
  range, th.ngals
  for ih = 0L, N_ELEMENTS(h) - 1 do begin
    if (h[i].ngals gt 0 and th[i].ngals gt 0) then h[i].edge = 1
    h[i].lbcg = min([h[i].lbcg, th[i].lbcg])
  endfor
  h.ngals += th.ngals
  h.n18 += th.n18
  h.n19 += th.n19
  h.n20 += th.n20
  h.n21 += th.n21
  h.n22 += th.n22
  h.lumtot += th.lumtot
  h.lum20 += th.lum20
endfor

ii = where(h.z gt zmin and h.z le zmax and h.mvir ge mmin)

mwrfits, h[ii], path+outname, /create

end
