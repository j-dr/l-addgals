pro split_halo_catalog_healpix, infile, outbase, nside

h = mrdfits(infile,1)

theta = (90.-h.dec)*!PI/180.
phi = h.ra*!PI/180.

ang2pix_ring, nside, theta, phi, ip
hist = histogram(ip, min=0, bin = 1, reverse_indices=ri)
for i = 0L, N_ELEMENTS(hist) do begin
  if (hist[i] le 50) then continue
  ind = ri[ri[i]+lindgen(hist[i])]
  
  outfile = outbase+'.'+strcompress(string(i))+'.fit'
  mwrfits, h[ind], outfile, /create
endfor

end

