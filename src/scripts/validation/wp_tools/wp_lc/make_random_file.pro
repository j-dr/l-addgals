
; number of points to generate over the volume of the full sphere
npoints = 10000000

; where we'll save everything
outdir = '/nfs/slac/g/ki/ki22/cosmo/jderose/addgals/catalogs/Fox/Catalog_v1.0/random/'
fbase = 'ranfile_'

; the parameters of our different volume limited samples
cspeed = 2.99792e5
mrmin = ['18', '19', '20', '21']
zmax = [12500.0/cspeed, 19250.0/cspeed, 31900.0/cspeed, 47650./cspeed]

; the geometery of our catalog
nside = 4
indices = ['0', '1', '4', '5', '6', '7']

; calibrate our redshift-distance relation
omegam = 0.31834
omegal = 0.68166
table = generate_z_of_r_table(omegam, omegal)

; generate our random positions
phi = randomu(seed, npoints)*2*!PI
theta = acos(randomu(seed, npoints)*2-1)
ang2pix_ring, nside, theta, phi, ip

; downsample to ones that are just in our region of interest
hist = histogram(ip, min = 0, reverse_indices = ri)
ikeep = -1
nkeep = 0
for i = 0L, N_ELEMENTS(indices) - 1 do begin
  if (hist[indices[i]] eq 0) then continue
  ind = ri[ri[indices[i]]+lindgen(hist[indices[i]])]
  ikeep = [ikeep, ind]
  nkeep += hist[indices[i]]
endfor
ikeep = ikeep[1:N_ELEMENTS(ikeep)-1]
phi = phi[ikeep]
theta = theta[ikeep]
print, "Number of random points: ", nkeep

; convert to ra and dec
ra = phi*180./!PI
dec = 90.0-(theta*180./!PI)

; loop through our different samples
for i = 0, N_ELEMENTS(mrmin) - 1 do begin

  ; setup our outfile
  fout = outdir+'/'+fbase+mrmin[i]+'.txt'
  
  ; get our max radius for this sample
  rmax = r_of_z(zmax[i], table)
  rmax3 = rmax^3

  ; generate a bunch of random radii and convert to redshifts
  rran = rmax*randomu(seed, nkeep)^(1./3)
 
  ; convert radii to resfhits
  czran = cspeed*z_of_r(rran, table)

  ; save
  openw,1,fout
  niceprintf,1,ra, dec, czran
  close, 1

endfor

end
