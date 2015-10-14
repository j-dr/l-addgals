pro make_full_bcc_halo_catalog, path, outname1, outname2, outname3, outbase, pixnum=pixnum

if not KEYWORD_SET(pixnum) then pixnum = ['0', '1', '4', '5', '6', '7', '12', '13', '14', '15', '16', '17', '24', '25', '26', '27', '28', '29', '30', '31', '40', '41', '42', '43', '44', '45', '46', '47', '48', '56', '57', '58', '59', '60', '61', '62', '63', '72', '73', '74', '75', '76', '77', '78', '79', '80', '88', '89', '90', '91', '92', '93', '94', '95']

;;;combine the individual box files internally
zmin = 0.0
zmax = 0.34
mmin = 3e12
make_bcc_halo_catalog, path, outname1, outname1, pixnum, zmin=zmin, zmax=zmax, mmin=mmin

zmin = 0.34
zmax = 0.9
mmin = 3e12
make_bcc_halo_catalog, path, outname2, outname2, pixnum, zmin=zmin, zmax=zmax, mmin=mmin
 
zmin = 0.9
zmax = 2.0
mmin = 2.4e13
make_bcc_halo_catalog, path, outname3, outname3, pixnum, zmin=zmin, zmax=zmax, mmin=mmin

;;;create one master catalog
f1=path+'/'+outname1
f2=path+'/'+outname2
f3=path+'/'+outname3
combine_bcc_halo_catalog, f1, f2, f3, outbase+'.fit'

;;;rotate and split the master catalog
h = mrdfits(outbase+'.fit',1)

ra = h.ra
dec = h.dec
dr = sqrt(h.halopx*h.halopx + h.halopy*h.halopy + h.halopz*h.halopz)
add_rotation_angle, ra,dec, ra_new, dec_new, rpx=rpx, rpy=rpy, rpz=rpz, dr=dr
h.ra = ra_new
h.dec = dec_new
h.halopx = rpx
h.halopy = rpy
h.halopz = rpz

;;;split the catalog
nside = 2
theta = (90.-dec)*!PI/180.
phi = ra*!PI/180.

ang2pix_ring, nside, theta, phi, ip
hist = histogram(ip, min=0, bin = 1, reverse_indices=ri)
for i = 0L, N_ELEMENTS(hist)-1 do begin
  if (hist[i] le 50) then continue
  ind = ri[ri[i]+lindgen(hist[i])]

  outfile = outbase+'_rotated.'+strcompress(string(i),/remove_all)+'.fit'
  mwrfits, h[ind], outfile, /create
endfor

print, 'Created the final Halo Catalog.'

end

