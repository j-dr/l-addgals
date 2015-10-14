pro repackage_qso;, infile, outfile

infile = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/QSO/Brazil/v1.0/randpos.cat'
outpath = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/QSO/Brazil/v1.0/'
outbase = 'ADDQSO_v1.0'

spawn, 'mkdir '+outpath+'/truth/'
spawn, 'mkdir '+outpath+'/obs/'

nside = 8

print, 'reading data...'
fmt = 'l,f,f,f,f,f,f,f,f,f'
readcol, infile, id, ra, dec, z, amag, mg, mr, mi, mz, my, format=fmt
nqso = 2*N_ELEMENTS(id)
nqso1 = N_ELEMENTS(id)
ind1 = lindgen(nqso1)
ind2 = ind1+nqso1

print, 'Transfering to truth structure....'
truth1 = get_bcc_truth_structure()
truth = replicate(truth1, nqso)
truth[ind1].id = id
truth[ind1].tra = ra
truth[ind1].tdec = dec
truth[ind1].ra = ra
truth[ind1].dec = dec
truth[ind1].z = z
truth[ind1].amag[2] = amag
truth[ind1].tmag[0] = mg
truth[ind1].tmag[1] = mr
truth[ind1].tmag[2] = mi
truth[ind1].tmag[3] = mz
truth[ind1].tmag[4] = my
truth[ind1].pqso = 1.0
truth[ind2].id = id+max(id)+1
truth[ind2].tra = ra+90.
truth[ind2].tdec = dec
truth[ind2].ra = ra+90.
truth[ind2].dec = dec
truth[ind2].z = z
truth[ind2].amag[2] = amag
truth[ind2].tmag[0] = mg
truth[ind2].tmag[1] = mr
truth[ind2].tmag[2] = mi
truth[ind2].tmag[3] = mz
truth[ind2].tmag[4] = my
truth[ind2].pqso = 1.0

;id = 0
;z = 0
;amag = 0
;mg = 0
;mr = 0
;mi = 0
;mz = 0
;my = 0

;;; add photometric errors
print, 'Adding Photometric Errors....'
mock_error_apply, 'DES', truth.omag, flux, ivar, omag, omagerr, /point_source
truth.omag = omag
truth.omagerr = omagerr
truth.flux = flux
truth.ivar = 1./ivar^2

;;;get our pixel decomposition
print, 'Getting healpix decomposition....'
theta = (90-truth.dec)*!PI/180.
phi = truth.ra*!PI/180.
ang2pix_ring, nside, theta, phi, ip
hist = histogram(ip, min = 0, bin = 1, reverse_indices = ri)

;;;rotate the catalog
print, 'doing rotation....'
add_rotation_angle, truth.ra, truth.dec, new_ra, new_dec
truth.ra = new_ra
truth.dec = new_dec
truth.tra = new_ra
truth.tdec = new_dec

;;;make the obs catalog
print, 'Making obs catalog....'
obs1 = get_bcc_obs_structure()
obs = replicate(obs1, nqso)
obs.id = truth.id
obs.ra = truth.ra
obs.dec = truth.dec
obs.mag_g = truth.omag[0]
obs.mag_r = truth.omag[1]
obs.mag_i = truth.omag[2]
obs.mag_z = truth.omag[3]
obs.mag_y = truth.omag[4]
obs.magerr_g = truth.omagerr[0]
obs.magerr_r = truth.omagerr[1]
obs.magerr_i = truth.omagerr[2]
obs.magerr_z = truth.omagerr[3]
obs.magerr_y = truth.omagerr[4]
obs.flux_g = truth.flux[0]
obs.flux_r = truth.flux[1]
obs.flux_i = truth.flux[2]
obs.flux_z = truth.flux[3]
obs.flux_y = truth.flux[4]
obs.ivar_g = truth.ivar[0]
obs.ivar_r = truth.ivar[1]
obs.ivar_i = truth.ivar[2]
obs.ivar_z = truth.ivar[3]
obs.ivar_y = truth.ivar[4]
obs.pqso = truth.pqso

;;;loop through the pixels, writing our output catalog
print, 'Saving the catalog....'
for i = 0L, N_ELEMENTS(hist)-1 do begin
  if (hist[i] eq 0) then continue
  ind = ri[ri[i]+lindgen(hist[i])]
  spix = strcompress(string(i),/remove_all)
  print, "Writing pixel "+spix
  toutfile = outpath+'/truth/'+outbase+'_truth.'+spix+'.fit'
  ooutfile = outpath+'/obs/'+outbase+'.'+spix+'.fit'
  mwrfits, truth[ind], toutfile, /create
  mwrfits, obs[ind], ooutfile, /create
endfor

end
