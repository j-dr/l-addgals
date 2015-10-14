;pro convert_cooper_data

; used for converting angular size to a physical size
omegam = 0.3
omegal = 0.7

;;;new g-r cut = 0.1 - 0.035*gmr.  Need to re-calculate the fraction

g = mrdfits('/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/cooper/sdss.environ.dr6.fits',1)
spec = mrdfits('/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/object_sdss_spectro.fits',1)
kgals = mrdfits('/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/kcorrect.photoz.model.z0.10.fits',1)
k = mrdfits('/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/object_sdss_imaging.fits',1)
sersic = mrdfits('/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/sersic_catalog.fits',1)
outfile = '/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/cooper/combined_dr6_cooper_extended.fit'
outtable = '/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/cooper/dr6_cooper_id_with_red_extended.dat'

zmax = 0.2
zmin = 0.0

zbin_size = 0.05
nzbins = floor((zmax-zmin)/zbin_size)
ninbin = lonarr(nzbins)

ii = where(g.z lt zmax)
g = g[ii]
ng = N_ELEMENTS(g)
percent5 = fltarr(ng)
percent6 = fltarr(ng)

for i = 0, nzbins - 1 do begin
  this_zmin = zmin + i*zbin_size
  this_zmax = this_zmin + zbin_size
  inbin = where(g.z gt this_zmin and g.z le this_zmax, tninbin)
  ninbin(i) = tninbin
  srt = sort(g(inbin).Delta5)
  for ii = 0L, tninbin-1 do begin
    percent5(inbin(srt(ii))) = (ii+0.5)/float(tninbin)
  endfor
  srt = sort(g(inbin).Delta6)
  for ii = 0L, tninbin-1 do begin
    percent6(inbin(srt(ii))) = (ii+0.5)/float(tninbin)
  endfor
  srt = 0
endfor

;super1 = g.plate+g.mjd*10LL^3+g.fiberid*10LL^7
;super2 = spec.plate+spec.mjd*10LL^3+spec.fiberid*10LL^7

in = where(spec.plate ge 0 and spec.z lt zmax and spec.z ge 0)
spec = spec[in]
k = k[in]
kgals = kgals[in]
sersic = sersic[in]
table = generate_z_of_r_table(omegam, omegal)
rlos = r_of_z(spec.z, table)
da = rlos/(1+spec.z)
theta = sersic.sersic_r0[2]/60/60*!PI/180.0
gsize = theta*da
;in = 0

plate_min = min([g.plate, spec.plate])
mjd_min = min([g.mjd, spec.mjd])
fiberid_min = min([g.fiberid, spec.fiberid])
g.plate -= plate_min
g.mjd -= mjd_min
g.fiberid -= fiberid_min
spec.plate -= plate_min
spec.mjd -= mjd_min
spec.fiberid -= fiberid_min
plate_max = long64(max([g.plate, spec.plate])+1)
mjd_max = long64(max([g.mjd, spec.mjd])+1)
fiberid_max = long64(max([g.fiberid, spec.fiberid])+1)
super1 = g.plate+g.mjd*plate_max+g.fiberid*plate_max*mjd_max
super2 = spec.plate+spec.mjd*plate_max+spec.fiberid*plate_max*mjd_max

sort1 = sort(super1)
sort2 = sort(super2)
i2 = 0L
matches = lonarr(ng)
for i1 = 0L, ng - 1 do begin
  while (super1(sort1(i1)) gt super2(sort2(i2))) do begin
    i2++
  endwhile
  if (super1(sort1(i1)) ne super2(sort2(i2))) then begin
    print, "ERROR!  super1 != super 2"
    stop
  endif
  matches(sort1(i1)) = sort2(i2)
endfor

;super1 = g.plate+g.mjd*10LL^3+g.fiberid*10LL^6
;super2 = spec.plate+spec.mjd*10LL^3+spec.fiberid*10LL^6
;match, super1, super2, a, b


out_arr = create_struct('index', 0L, 'amag',fltarr(5), 'coeffs',fltarr(5), 'Delta5', 0., 'Delta6', 0.,$
                        'PercentDensity5', 0., 'PercentDensity6', 0., 'kc', fltarr(5), 'z',0.,$
                       'abmaggies', fltarr(5), 'sdss_environ_id', 0L, $
			'mass', 0., 'intsfh', 0., 'mets', 0., 'b300', 0., 'b1000', 0., $
			'n_sersic', 0., 'r0_sersic', 0.)
gals = replicate(out_arr, ng)



kcorrect = sdss_kcorrect(spec.z, calibobj=k, absmag=absmag, band_shift=0.1, $
                         flux="model", coeffs=coeffs, nmgy=nmgy, ivar=ivar)
this_mag = fltarr(ng)
this_coeffs = fltarr(ng)
gals.index = in[matches]
for i = 0, 4 do begin
;  gals(*).amag(i) = k(matches).absmag(i)
;  gals(*).amag(i) = g(*).absmag(i)
  this_mag(*) = absmag(i,matches)
  gals(*).amag(i) = this_mag(*)
  this_coeffs(*) = coeffs(i,matches)
  gals(*).coeffs(i) = this_coeffs(*)
  gals(*).abmaggies(i) = kgals(matches).abmaggies(i)
endfor
gals.n_sersic = sersic[matches].sersic_n[2]
gals.r0_sersic = sersic[matches].sersic_r50[2]
gals.Delta5 = g.Delta5
gals.Delta6 = g.Delta6
gals.PercentDensity5 = percent5
gals.PercentDensity6 = percent6
gals.kc = g.kcorrect
gals.z = g.z
;gals(*).coeffs = k(matches).coeffs
gals.sdss_environ_id = matches
gals.mass = kgals[matches].mass
gals.intsfh = kgals[matches].intsfh
gals.mets = kgals[matches].mets
gals.b300 = kgals[matches].b300
gals.b1000 = kgals[matches].b1000

gmr = gals.amag(1) - gals.amag(2)
ii = where(gmr le 3.0)
gals = gals[ii]

stop

mwrfits, gals, outfile, /create

gmr = gals.amag(1) - gals.amag(2)
red = intarr(N_ELEMENTS(gals))
red(*) = 0
;ii = where(gmr ge 0.72)
ii = where (gmr ge 0.1 - 0.035*gals.amag[2])
red(ii) = 1

openw,1,outtable
for i = 0L, N_ELEMENTS(gals)-1 do begin
  printf,1,gals(i).amag(2), gals(i).PercentDensity6, gals(i).PercentDensity5, red(i)
endfor
close,1

end
