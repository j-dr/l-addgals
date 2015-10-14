pro combine_healpix_catalogs, inpath, fbase, outpath, sizes, nside_in, nside_out

nfiles = N_ELEMENTS(sizes)
npix = 12*nside_in*nside_in

truth_base = inpath+fbase+'_'+sizes+'_truth_no_photoz.'
obs_base = inpath+fbase+'_'+sizes+'_no_photoz.'
deep_base = inpath+fbase+'_'+sizes+'_deep_mag.'
johnson_base = inpath+fbase+'_'+sizes+'_johnson_mag.'
flamex_base = inpath+fbase+'_'+sizes+'_flamex_mag.'
cfhtls_base = inpath+fbase+'_'+sizes+'_cfhtls_mag.'
vista_base = inpath+fbase+'_'+sizes+'_vista_mag.'
sdss_base = inpath+fbase+'_'+sizes+'_sdss_mag.'

otruth_base = outpath+fbase+'_truth_no_photoz.'
oobs_base = outpath+fbase+'_no_photoz.'
odeep_base = outpath+fbase+'_deep_mag.'
ojohnson_base = outpath+fbase+'_johnson_mag.'
oflamex_base = outpath+fbase+'_flamex_mag.'
ocfhtls_base = outpath+fbase+'_cfhtls_mag.'
ovista_base = outpath+fbase+'_vista_mag.'
osdss_base = outpath+fbase+'_sdss_mag.'


;for i = 0L, npix - 1 do begin
for i = 4, 4 do begin
  sfx = strcompress(string(i), /remove_all)+'.fit'
  intruth = truth_base+sfx
  inobs = obs_base+sfx
  indeep = deep_base+sfx
  injohnson = johnson_base+sfx
  inflamex = flamex_base+sfx
  invista = vista_base+sfx
  insdss = sdss_base+sfx
  if (file_test(intruth[0]) eq 0) then continue 
  for j = 0, nfiles - 1 do begin
    truth = mrdfits(intruth[j],1)
    obs = mrdfits(inobs[j],1)
    deep = mrdfits(indeep[j],1)
    johnson = mrdfits(injohnson[j],1)
    flamex=mrdfits(inflamex[j],1)
    vista=mrdfits(invista[j],1)
    sdss=mrdfits(insdss[j],1)
    
    theta = (90-truth.dec)*!PI/180.
    phi = truth.ra*!PI/180.
    ang2pix_ring, nside_out, theta, phi, ip_out 

    hist = histogram(ip_out, min = 0, reverse_indices = ri)
    for ih = 0L, N_ELEMENTS(hist) - 1 do begin
      sfx2 = strcompress(string(ih), /remove_all)+'.fit'
      if (hist[ih] lt 10) then continue
      ii = ri[ri[ih]:ri[ih]+hist[ih]-1]
      outtruth = otruth_base+sfx2
      outobs = oobs_base+sfx2
      outdeep = odeep_base+sfx2
      outjohnson = ojohnson_base+sfx2
      outflamex = oflamex_base+sfx2
      outvista = ovista_base+sfx2
      outsdss = osdss_base+sfx2

;      append = (FILE_TEST(outtruth) eq 1) ? 1 : 0
;      create = (append eq 0) ? 1 : 0
      
;      create = (j eq 0) ? 1 : 0

;      mwrfits, truth[ii], outtruth, create=create
;      mwrfits, obs[ii], outobs, create=create
;      mwrfits, deep[ii], outdeep,create=create
;      mwrfits, johnson[ii], outjohnson, create=create
;      mwrfits, flamex[ii], outflamex, create=create
;      mwrfits, vista[ii], outvista, create=create
;      mwrfits, sdss[ii], outsdss, create=create

;      if (0) then begin

      if (j eq 0) then begin
        mwrfits, truth[ii], outtruth, /create
        mwrfits, obs[ii], outobs, /create
        mwrfits, deep[ii], outdeep, /create
        mwrfits, johnson[ii], outjohnson, /create
        mwrfits, flamex[ii], outflamex, /create
        mwrfits, vista[ii], outvista, /create
        mwrfits, sdss[ii], outsdss, /create
      endif else begin
        ttruth = mrdfits(outtruth,1)
	tobs = mrdfits(outobs,1)
	tdeep = mrdfits(outdeep,1)
	tjohnson = mrdfits(outjohnson,1)
	tflamex = mrdfits(outflamex,1)
	tvista = mrdfits(outvista,1)
	tsdss = mrdfits(outsdss,1)

        mwrfits, [ttruth, truth[ii]], outtruth, /create
        mwrfits, [tobs, obs[ii]], outobs, /create
        mwrfits, [tdeep, deep[ii]], outdeep, /create
        mwrfits, [tjohnson, johnson[ii]], outjohnson, /create
        mwrfits, [tflamex, flamex[ii]], outflamex, /create
        mwrfits, [tvista, vista[ii]], outvista, /create
        mwrfits, [tsdss, sdss[ii]], outsdss, /create
      endelse

;      endif
   
    endfor
  endfor
endfor

end

