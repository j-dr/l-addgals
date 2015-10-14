pro reformat_addstars_file, infile, outpath, nside=nside

strarr = strsplit(infile, '/', /extract)
fname = strarr[N_ELEMENTS(strarr)-1]

s = mrdfits(infile,1)
obs1 = get_bcc_obs_structure()
truth1 = get_bcc_truth_structure()
remove_rotation_angle, s.ra, s.dec, ra, dec

if not KEYWORD_SET(nside) then nside = 8
npix = 12*nside*nside
theta = (90. - dec)*!PI/180.
phi = ra*!PI/180.

ang2pix_ring, nside, theta, phi, ip
hist = histogram(ip, min = 0, reverse_indices = ri)

for ih = 0L, N_ELEMENTS(hist) - 1 do begin
  if(hist[ih] eq 0) then continue
  ii = ri[ri[ih]:ri[ih]+hist[ih]-1]
  num = strcompress(string(ih),/remove_all)
;  truth_file = outpath+'/PO_Aardvark_stars.'+num+'.fit'
;  obs_file = outpath+'/PO_Aardvark_truth_stars.'+num+'.fit'
  truth_file = outpath+'/'+fname+'_truth.'+num+'.fit'
  obs_file = outpath+'/'+fname+'_obs.'+num+'.fit'

  obs = replicate(obs1,hist[ih])
  truth = replicate(truth1,hist[ih])

;  truth.id = 
  truth.tmag[0] = s[ii].gdes
  truth.tmag[1] = s[ii].rdes
  truth.tmag[2] = s[ii].ides
  truth.tmag[3] = s[ii].zdes
  truth.tmag[4] = s[ii].ydes
  mock_error_apply,'DES',truth.tmag,flux,ivar,omag,omagerr
  truth.omag = omag
  truth.flux = flux
  truth.ivar = 1./ivar^2
  truth.omagerr = omagerr
  truth.ra = s[ii].ra
  truth.dec = s[ii].dec
  truth.haloid = -1
  truth.tra = s[ii].ra
  truth.tdec = s[ii].dec
  truth.lmag = truth.tmag
  truth.mag_u = s[ii].usds
  truth.pstar = 1.0

  obs.id = truth.id
  obs.mag_u = s[ii].usds
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
  obs.ra = truth.ra
  obs.dec = truth.dec
  obs.pstar = 1.0
  
;  if (file_test(truth_file)) then begin
;    in_truth = mrdfits(truth_file, 1)
;    truth = [in_truth, truth]
;    in_obs = mrdfits(obs_file,1)
;    obs = [in_obs, obs]
;  endif
  mwrfits, truth, truth_file, /create
  mwrfits, obs, obs_file, /create
endfor


end
