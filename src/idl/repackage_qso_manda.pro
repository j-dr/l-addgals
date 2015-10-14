pro repackage_qso_manda;, infile, outfile

infile = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/DESQSO/DES_VHS_Allz_Mock_20120521_mbanerji'
outpath = '/nfs/slac/g/ki/ki19/des/mbusha/catalogs/QSO/DESQSO/'
outbase = 'DESQSO'
add_lensing = 0

spawn, 'mkdir '+outpath+'/truth/'
spawn, 'mkdir '+outpath+'/obs/'
spawn, 'mkdir '+outpath+'/VHS/'
;spawn, 'mkdir '+outpath+'/VIKING/'

nside = 8

print, 'reading data...'
fmt = 'l,f,f,f,f,f,f,f,f,f'
rdfloat, infile, z, mg, mr, mi, mz, my, mj, mh, mk
nqso = N_ELEMENTS(z)
id = lindgen(nqso)

print, 'Transfering to truth structure....'
truth1 = get_bcc_truth_structure()
truth = replicate(truth1, nqso)
truth.id = id
truth.tra = 180.*randomu(seed, nqso)
truth.tdec = asin(randomu(seed, nqso))*180./!PI
truth.ra = truth.tra
truth.dec = truth.tdec
truth.z = z
truth.tmag[0] = mg
truth.tmag[1] = mr
truth.tmag[2] = mi
truth.tmag[3] = mz
truth.tmag[4] = my
truth.lmag = truth.tmag
truth.pqso = 1.0

;;; add photometric errors
print, 'Adding Photometric Errors....'
mock_error_apply, 'DES', truth.omag, flux, ivar, omag, omagerr, /point_source
truth.omag = omag
truth.omagerr = omagerr
truth.flux = flux
truth.ivar = 1./ivar^2

;;;make the VHS/VIKING catalogs
vhs1 = create_struct('tmag', fltarr(3),'lmag',fltarr(3),'amag',fltarr(3), $
                     'omag', fltarr(3), 'omagerr', fltarr(3), $
                     'flux',fltarr(3), 'ivar', fltarr(3))
vhs = replicate(vhs1,nqso)
vhs.tmag[0] = mj
vhs.tmag[1] = mh
vhs.tmag[2] = mk
vhs.lmag = vhs.tmag
print, "adding vhs photometric errors...."
if (add_lensing) then begin
  mock_error_apply,'VHS',vhs.lmag,flux,ivar,omag,omagerr
endif else begin
  mock_error_apply,'VHS',vhs.tmag,flux,ivar,omag,omagerr
endelse
vhs.omag = omag
vhs.omagerr = omagerr
vhs.flux = flux
vhs.ivar = 1./ivar^2

;print, "adding viking photometric errors...."
;viking = vhs
;if (add_lensing) then begin
;  mock_error_apply,'VIKING',viking.lmag,flux,ivar,omag,omagerr
;endif else begin
;  mock_error_apply,'VIKING',viking.tmag,flux,ivar,omag,omagerr
;endelse
;viking.omag = omag
;viking.omagerr = omagerr
;viking.flux = flux
;viking.ivar = 1./ivar^2

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
  vhsoutfile = outpath+'/VHS/'+outbase+'.'+spix+'.fit'
  vikingoutfile = outpath+'/VIKING/'+outbase+'.'+spix+'.fit'
  mwrfits, truth[ind], toutfile, /create
  mwrfits, obs[ind], ooutfile, /create
  mwrfits, vhs[ind], vhsoutfile, /create
;  mwrfits, viking[ind], vikingoutfile, /create
endfor

end
