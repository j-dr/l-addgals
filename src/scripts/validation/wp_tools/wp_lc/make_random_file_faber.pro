pro make_random_file_sva
; number of points to generate over the volume of the full sphere
npoints = 10000000

; where we'll save everything
outdir = '/nfs/slac/g/ki/ki23/des/jderose/addgals/catalogs/chinchilla-tuning/SVA_faber_lf_evol/validation/random'
fbase = 'ranfile_'

; the parameters of our different volume limited samples
cspeed = 2.99792e5
mrmin = ['18', '19', '20', '21']
zmax = [12500.0/cspeed, 19250.0/cspeed, 31900.0/cspeed, 47650./cspeed]

; the geometery of our catalog
nside = 32
indices = ['10074', '10905', '10683', '10448', '11202', '10206', '10907',$
 '10687', '10080', '10456', '11296', '10327', '10903', '10573', '10078',$
 '10452', '11200', '10333', '10801', '10575', '10200', '10454', '11378',$
 '10331', '10799', '10795', '10198', '10569', '11294', '10070', '10689',$
 '10909', '10329', '11012', '10571', '11206', '10072', '10691', '10901',$
 '11105', '10202', '10567', '11204', '10076', '11008', '10797', '10685',$
 '10325', '10450', '11006', '11292', '11107', '10204', '11010', '11111',$
 '11109', '11462', '11464', '11384', '11382', '11466', '11380', '10205',$
 '10207', '10326', '10330', '10332', '10203', '11203', '10334', '10328',$
 '10568', '10453', '11295', '11291', '11205', '10570', '10684', '11293',$
 '11379', '11297', '10451', '10457', '10449', '10576', '10572', '11201',$
 '10455', '10796', '10794', '10686', '10688', '10690', '10574', '10199',$
 '10075', '10197', '10079', '10201', '10073', '10071', '10077', '10798',$
 '10802', '11005', '10904', '10902', '10908', '10906', '10800', '11465',$
 '11461', '11467', '11383', '11381', '11463', '11106', '11009', '11011',$
 '11104', '11110', '11108', '11199', '11007']

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
  fout = outdir+'/'+fbase+mrmin[i]+'.dat'
  
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
