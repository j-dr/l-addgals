pro combine_single_healpix_catalog, pixel, inpath, shearbase, fbase, outpath, $
	sizes, nside_in, nside_out, mag_lim=mag_lim, obase=obase, $
	add_des_errors=add_des_errors, $
	add_vista_errors=add_vista_errors, $
	add_dr8_errors=add_dr8_errors, $
	add_stripe82_errors=add_stripe82_errors, $
	add_rcs_errors=add_rcs_errors, $
	add_cfhtls_errors=add_cfhtls_errors, $
	add_deep_errors=add_deep_errors, $
	add_flamex_errors=add_flamex_errors, $
	add_lensing=add_lensing

;;;depth to which we cut the catalog
sigma_cut = 5.0

;;;set some healpix definitions
print, "combining pixel ", pixel
nmingals = 1000
nfiles = N_ELEMENTS(sizes)
npix = 12*nside_in*nside_in
if not KEYWORD_SET(obase) then obase = fbase

boxinpath = inpath+'/individual_box_files/'

;;;define out input file base (add sfx for healpix number later)
truth_base = inpath+'/shapes/'+fbase+'_'+sizes+'_truth_no_photoz.'
obs_base = boxinpath+fbase+'_'+sizes+'_no_photoz.'
deep_base = boxinpath+fbase+'_'+sizes+'_deep_mag.'
johnson_base = boxinpath+fbase+'_'+sizes+'_johnson_mag.'
flamex_base = boxinpath+fbase+'_'+sizes+'_flamex_mag.'
vista_base = boxinpath+fbase+'_'+sizes+'_vista_mag.'
sdss_base = boxinpath+fbase+'_'+sizes+'_sdss_mag.'
cfhtls_base = boxinpath+fbase+'_'+sizes+'_cfhtls_mag.'
euclid_base = boxinpath+fbase+'_'+sizes+'_euclid_mag.'
irac_base = boxinpath+fbase+'_'+sizes+'_irac_mag.'
wise_base = boxinpath+fbase+'_'+sizes+'_wise_mag.'
hsc_base = boxinpath+fbase+'_'+sizes+'_hsc_mag.'
lsst_base = boxinpath+fbase+'_'+sizes+'_lsst_mag.'
wfirst_base = boxinpath+fbase+'_'+sizes+'_wfirst_mag.'
twomass_base = boxinpath+fbase+'_'+sizes+'_twomass_mag.'
;shear_base = shearpath+'/shearcat_'+fbase+'_'+sizes+'_truth_no_photoz.'
;shear_base = shearpath+'/shearcat_PO_Aardvark_'+sizes+'_truth_no_photoz.'
shear_base = shearbase+'_'+sizes+'_parts_'

;;;make out output directories
spawn, 'mkdir -p '+outpath+'/truth/'
spawn, 'mkdir -p '+outpath+'/obs/'
spawn, 'mkdir -p '+outpath+'/DEEP2/'
spawn, 'mkdir -p '+outpath+'/Johnson/'
spawn, 'mkdir -p '+outpath+'/FLAMEX/'
spawn, 'mkdir -p '+outpath+'/NDWFS/'
spawn, 'mkdir -p '+outpath+'/RCS/'
spawn, 'mkdir -p '+outpath+'/WISE/'
spawn, 'mkdir -p '+outpath+'/CFHTLS/'
spawn, 'mkdir -p '+outpath+'/VHS/'
spawn, 'mkdir -p '+outpath+'/VIKING/'
spawn, 'mkdir -p '+outpath+'/DR8/'
spawn, 'mkdir -p '+outpath+'/Stripe82/'
spawn, 'mkdir -p '+outpath+'/Euclid/'
spawn, 'mkdir -p '+outpath+'/HSC/'
spawn, 'mkdir -p '+outpath+'/LSST/'
spawn, 'mkdir -p '+outpath+'/WFIRST/'
spawn, 'mkdir -p '+outpath+'/2MASS/'
spawn, 'mkdir -p '+outpath+'/BCS/'
spawn, 'mkdir -p '+outpath+'/IRAC/'
spawn, 'mkdir -p '+outpath+'/WISE/'

;;;define the output files
otruth_base = outpath+'/truth/'+obase+'_truth.'
oobs_base = outpath+'/obs/'+obase+'.'
odeep_base = outpath+'/DEEP2/'+obase+'_deep_mag.'
ojohnson_base = outpath+'/Johnson/'+obase+'_johnson_mag.'
oflamex_base = outpath+'/FLAMEX/'+obase+'_flamex_mag.'
ondwfs_base = outpath+'/NDWFS/'+obase+'_ndwfs_mag.'
orcs_base = outpath+'/RCS/'+obase+'_rcs_mag.'
owise_base = outpath+'/WISE/'+obase+'_wise_mag.'
ocfhtls_base = outpath+'/CFHTLS/'+obase+'_cfhtls_mag.'
ovista_base = outpath+'/VHS/'+obase+'_vista_mag.'
ovhs_base = outpath+'/VHS/'+obase+'_vhs_mag.'
oviking_base = outpath+'/VIKING/'+obase+'_viking_mag.'
oeuclid_base = outpath+'/Euclid/'+obase+'_euclid_mag.'
oirac_base = outpath+'/IRAC/'+obase+'_irac_mag.'
owise_base = outpath+'/WISE/'+obase+'_wise_mag.'
ohsc_base = outpath+'/HSC/'+obase+'_hsc_mag.'
olsst_base = outpath+'/LSST/'+obase+'_lsst_mag.'
owfirst_base = outpath+'/WFIRST/'+obase+'_wfirst_mag.'
otwomass_base = outpath+'/2MASS/'+obase+'_2mass_mag.'
osdss_base = outpath+'/DR8/'+obase+'_sdss_mag.'
ostripe82_base = outpath+'/Stripe82/'+obase+'_stripe82_mag.'
obcs_base = outpath+'/BCS/'+obase+'_bcs_mag.'

  ;;;add the pixel information to definie explicit input files
  i = pixel
  sfx = i+'.fit'
  intruth = truth_base+sfx
  inobs = obs_base+sfx
  indeep = deep_base+sfx
  injohnson = johnson_base+sfx
  inflamex = flamex_base+sfx
  incfhtls = cfhtls_base+sfx
  invista = vista_base+sfx
  ineuclid = euclid_base+sfx
  inirac = irac_base+sfx
  inwise = wise_base+sfx
  inhsc = hsc_base+sfx
  inlsst = lsst_base+sfx
  inwfirst = wfirst_base+sfx
  intwomass = twomass_base+sfx
  insdss = sdss_base+sfx
  inshear = shear_base+sfx

  ;;;loop through the different box files
  for j = 0, nfiles - 1 do begin
    print, "Reading in files..."
    print, intruth[j]
    print, inshear[j]

    ;;;always read the truth file
    truth = mrdfits(intruth[j],1)

    ;;;read the shear catalog for lensing and lens positions, if specified
    if KEYWORD_SET(add_lensing) then begin
      shear = mrdfits(inshear[j],1)
      print, 'adding lensing to truth catalog....'
      add_lensing, truth, shear
    endif

    ;;;add photometric errors for DES catalog system
    if KEYWORD_SET(add_des_errors) then begin
      print, 'adding DES photometric errors...'
      if (add_lensing) then begin
        help, /struct, truth
        mock_error_apply,'DES',truth.lmag,flux,ivar,omag,omagerr
      endif else begin
        mock_error_apply,'DES',truth.tmag,flux,ivar,omag,omagerr
      endelse
      truth.omag = omag
      truth.omagerr = omagerr
      truth.flux = flux
      truth.ivar = 1./ivar^2
    endif

    ;;;cut at the 10-sigma limit
    print, 'cutting to 10-sigma'
    ind = where(truth.flux[0]*sqrt(truth.ivar[0]) gt sigma_cut or $
		   truth.flux[1]*sqrt(truth.ivar[1]) gt sigma_cut or $
		   truth.flux[2]*sqrt(truth.ivar[2]) gt sigma_cut or $
		   truth.flux[3]*sqrt(truth.ivar[3]) gt sigma_cut or $
		   truth.flux[4]*sqrt(truth.ivar[4]) gt sigma_cut, ng)
    print, "Fraction of galaxies to be kept: ", float(N_ELEMENTS(ind))/float(N_ELEMENTS(truth))
;    srt = sort(truth[ind].id)
;    u = uniq(truth[ind].id, srt)
;    ind = ind[u]
    truth = truth[ind] 

    ;;;read the deep catalog, lens, and add photometric errors.  
    deep = mrdfits(indeep[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, deep, shear, 3
    deep = deep[ind]
    print, "adding deep photometric errors...."
    if (add_lensing) then begin
      mock_error_apply,'DEEP2',deep.lmag,flux,ivar,omag,omagerr
    endif else begin
      mock_error_apply,'DEEP2',deep.tmag,flux,ivar,omag,omagerr
    endelse
    add_tags, deep, ['omag', 'omagerr', 'flux', 'ivar'], ['fltarr(3)', 'fltarr(3)', 'fltarr(3)', 'fltarr(3)'], deep2
    deep2.omag = omag
    deep2.omagerr = omagerr
    deep2.flux = flux
    deep2.ivar = 1./ivar^2
    deep = deep2


    ;;;readin the deep catalog and lens -- no photometric errors for this system
    johnson = mrdfits(injohnson[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, johnson, shear, 5
    johnson = johnson[ind]

    ;;;readin the CFHTLS catalog, lens, then make CFHTLS and RCS versions.
    cfhtls=mrdfits(incfhtls[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, cfhtls, shear, 5
    cfhtls = cfhtls[ind]

    print, "adding cfhtls photometric errors...."
    if (add_lensing) then begin
      mock_error_apply,'CFHTLS',cfhtls.lmag,flux,ivar,omag,omagerr
    endif else begin
      mock_error_apply,'CFHTLS',cfhtls.tmag,flux,ivar,omag,omagerr
    endelse
    add_tags, cfhtls, ['omag', 'omagerr', 'flux', 'ivar'], ['fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], cfhtls2
    cfhtls2.omag = omag
    cfhtls2.omagerr = omagerr
    cfhtls2.flux = flux
    cfhtls2.ivar = 1./ivar^2
    cfhtls = cfhtls2

    print, "adding rcs photometric errors...."
    if (add_lensing) then begin
      mock_error_apply,'RCS',cfhtls.lmag[1:4],flux,ivar,omag,omagerr
    endif else begin
      mock_error_apply,'RCS',cfhtls.tmag[1:4],flux,ivar,omag,omagerr
    endelse
    rcs1 = create_struct('tmag', fltarr(4), 'lmag', fltarr(4), 'amag', fltarr(4), $
                          'omag', fltarr(4), 'omagerr', fltarr(4), $
                          'flux',fltarr(4), 'ivar', fltarr(4))
    rcs = replicate(rcs1,ng)
    rcs.tmag = cfhtls.tmag[1:4]
    if (tag_exist(cfhtls, 'lmag')) then rcs.lmag = cfhtls.lmag[1:4]
    rcs.amag = cfhtls.amag[1:4]
    rcs.omag = omag
    rcs.omagerr = omagerr
    rcs.flux = flux
    rcs.ivar = 1./ivar^2

    ;;;readin vista catalog, make both VHS and VIKING catalogs
    vista=mrdfits(invista[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, vista, shear, 5
    vista = vista[ind]

    print, "adding vhs photometric errors...."
    if (add_lensing) then begin
      mock_error_apply,'VHS',vista.lmag[2:4],flux,ivar,omag,omagerr
    endif else begin
      mock_error_apply,'VHS',vista.tmag[2:4],flux,ivar,omag,omagerr
    endelse
    vhs1 = create_struct('tmag', fltarr(3),'lmag',fltarr(3),'amag',fltarr(3), $
                          'omag', fltarr(3), 'omagerr', fltarr(3), $
                          'flux',fltarr(3), 'ivar', fltarr(3))
    vhs = replicate(vhs1,ng)
    vhs.tmag = vista.tmag[2:4]
    if (tag_exist(vista, 'lmag')) then vhs.lmag = vista.lmag[2:4]
    vhs.amag = vista.amag[2:4]
    vhs.omag = omag
    vhs.omagerr = omagerr
    vhs.flux = flux
    vhs.ivar = 1./ivar^2

    print, "adding viking photometric errors...."
    if (add_lensing) then begin
      mock_error_apply,'VIKING',vista.lmag,flux,ivar,omag,omagerr
    endif else begin
      mock_error_apply,'VIKING',vista.tmag,flux,ivar,omag,omagerr
    endelse
    add_tags, vista, ['omag', 'omagerr', 'flux', 'ivar'], ['fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], viking
    viking.omag = omag
    viking.omagerr = omagerr
    viking.flux = flux
    viking.ivar = 1./ivar^2
    vista = 0

    ;;;read in euclid, lens, but no photometric errors
    euclid=mrdfits(ineuclid[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, euclid, shear, 4
    euclid=euclid[ind]

    ;;;read in irac, lens, but no photometric errors
    irac=mrdfits(inirac[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, irac, shear, 4
    irac=irac[ind]

    ;;;read in wise, lens, but no photometric errors
    wise=mrdfits(inwise[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, wise, shear, 4
    wise=wise[ind]

    ;;;read in HSC, lens, no photometric errors
    hsc=mrdfits(inhsc[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, hsc, shear, 5
    hsc=hsc[ind]

    ;;;read in LSST, lens, no photometric errors
    lsst=mrdfits(inlsst[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, lsst, shear, 8
    lsst=lsst[ind]

    ;;;read in WFIRST, lens, no photometric errors
    wfirst=mrdfits(inwfirst[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, wfirst, shear, 4
    wfirst=wfirst[ind]

    ;;;twomass, lens, no photometric errors
    twomass=mrdfits(intwomass[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, twomass, shear, 3
    twomass=twomass[ind]


    sdss=mrdfits(insdss[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, sdss, shear, 5
    sdss = sdss[ind]

;    if KEYWORD_SET(add_dr8_errors) then begin
      print, "adding dr8 photometric errors...."
      if (add_lensing) then begin
	mock_error_apply,'DR8',sdss.lmag,flux,ivar,omag,omagerr
      endif else begin
	mock_error_apply,'DR8',sdss.tmag,flux,ivar,omag,omagerr
      endelse
      add_tags, sdss, ['omag', 'omagerr', 'flux', 'ivar'], ['fltarr(5)', 'fltarr(5)', 'fltarr(5)', 'fltarr(5)'], sdss2
      sdss2.omag = omag
      sdss2.omagerr = omagerr
      sdss2.flux = flux
      sdss2.ivar = 1./ivar^2
      sdss = sdss2
;    endif

;    if KEYWORD_SET(add_stripe82_errors) then begin
      print, "adding stripe82 photometric errors...."
      if (add_lensing) then begin
        mock_error_apply,'STRIPE82',sdss.lmag,flux,ivar,omag,omagerr
      endif else begin
        mock_error_apply,'STRIPE82',sdss.tmag,flux,ivar,omag,omagerr
      endelse
      stripe82 = sdss
      stripe82.omag = omag
      stripe82.omagerr = omagerr
      stripe82.flux = flux
      stripe82.ivar = 1./ivar^2
;    endif

    ;;;add bcs errors
    print, "Adding BCS photometric errors...."
    if (add_lensing) then begin
      mock_error_apply,'BCS',sdss.lmag[1:4],flux,ivar,omag,omagerr
    endif else begin
      mock_error_apply,'BCS',sdss.tmag[1:4],flux,ivar,omag,omagerr
    endelse

    bcs1 = create_struct('tmag', fltarr(4), 'lmag', fltarr(4), 'amag', fltarr(4), $
                         'omag', fltarr(4), 'omagerr', fltarr(4), $
                         'flux',fltarr(4), 'ivar', fltarr(4))
    bcs = replicate(bcs1,ng)
    bcs.amag = sdss.amag[1:4]
    bcs.tmag = sdss.tmag[1:4]
    if (tag_exist(sdss, 'lmag')) then bcs.lmag = sdss.lmag[1:4]
    bcs.omag = omag
    bcs.omagerr = omagerr
    bcs.flux = flux
    bcs.ivar = 1./ivar^2

    ;;;readin flamex catalog, lens, make NDWFS, FLAMEX, WISE catralogs
    flamex=mrdfits(inflamex[j],1)
    if KEYWORD_SET(add_lensing) then lens_magnitudes, flamex, shear, 8
    flamex = flamex[ind]

    print, "Adding flamex errors...."
    if (add_lensing) then begin
      mock_error_apply,'NDWFS',flamex.lmag[0:2],flux1,ivar1,omag1,omagerr1
      mock_error_apply,'FLAMEX',flamex.lmag[[3,5]],flux2,ivar2,omag2,omagerr2
      mock_error_apply,'WISE',flamex.lmag[6:7],flux3,ivar3,omag3,omagerr3
    endif else begin
      mock_error_apply,'NDWFS',flamex.tmag[0:2],flux1,ivar1,omag1,omagerr1
      mock_error_apply,'FLAMEX',flamex.tmag[[3,5]],flux2,ivar2,omag2,omagerr2
      mock_error_apply,'WISE',flamex.tmag[6:7],flux3,ivar3,omag3,omagerr3
    endelse

    ndwfs1 = create_struct('tmag', fltarr(3), 'lmag', fltarr(3), 'amag', fltarr(3), $
                           'omag', fltarr(3), 'omagerr', fltarr(3), $
                           'flux',fltarr(3), 'ivar', fltarr(3))
    ndwfs = replicate(ndwfs1,ng)
    ndwfs.amag = flamex.amag[0:2]
    ndwfs.tmag = flamex.tmag[0:2]
    if (tag_exist(flamex, 'lmag')) then ndwfs.lmag = flamex.lmag[0:2]
    ndwfs.omag = omag1
    ndwfs.omagerr = omagerr1
    ndwfs.flux = flux1
    ndwfs.ivar = 1./ivar1^2

    flamex1 = create_struct('tmag', fltarr(2), 'lmag', fltarr(2), 'amag', fltarr(2), $
                          'omag', fltarr(2), 'omagerr', fltarr(2), $
                          'flux',fltarr(2), 'ivar', fltarr(2))
    flamex_out = replicate(flamex1,ng)
    flamex_out.amag = flamex.amag[[3,5]]
    flamex_out.tmag = flamex.tmag[[3,5]]
    if (tag_exist(flamex, 'lmag')) then flamex_out.lmag = flamex.lmag[[3,5]]
    flamex_out.omag = omag2
    flamex_out.omagerr = omagerr2
    flamex_out.flux = flux2
    flamex_out.ivar = 1./ivar2^2

    wise1 = create_struct('tmag', fltarr(2), 'lmag', fltarr(2), 'amag', fltarr(2), $
                          'omag', fltarr(2), 'omagerr', fltarr(2), $
                          'flux',fltarr(2), 'ivar', fltarr(2))
    wise = replicate(wise1,ng)
    wise.amag = flamex.amag[6:7]
    wise.tmag = flamex.tmag[6:7]
    if (tag_exist(flamex, 'lmag')) then wise.lmag = flamex.lmag[6:7]
    wise.omag = omag3
    wise.omagerr= omagerr3
    wise.flux = flux3
    wise.ivar = 1./ivar3^2
    flamex = 0

    ;;;make sure galaxy id's are unique
    truth.id += j*1000000000LL

    ;;;re-make the obs file
    print, "creating observational catalog..."
    obs1 = get_bcc_obs_structure()
    obs = replicate(obs1, N_ELEMENTS(truth))
    obs.id = truth.id
    obs.index = truth.index
    obs.mag_u = truth.mag_u
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
    obs.epsilon1 = truth.epsilon[0]
    obs.epsilon2 = truth.epsilon[1]
    obs.size = truth.size
    obs.arborz = truth.arborz
    obs.arborz_err = truth.arborz_err
    obs.annz = truth.annz_err
    obs.photoz_gaussian = truth.photoz_gaussian

    ;;;calculate the new pixels for the repackaged catalog
    theta = (90-truth.tdec)*!PI/180.
    phi = truth.tra*!PI/180.
    ang2pix_ring, nside_out, theta, phi, ip_out 
    hist = histogram(ip_out, min = 0, reverse_indices = ri)

    ;;;loop through pixels and save all the catalogs
    print, hist
    for ih = 0L, N_ELEMENTS(hist) - 1 do begin
      if (hist[ih] lt nmingals) then continue
      sfx2 = strcompress(string(ih), /remove_all)+'.fit'
      print, "outputing ", hist[ih], " galaxies into pixel ", ih
      ;ii = ri[ri[ih]:ri[ih]+hist[ih]-1]
      ii = ri[ri[ih]+lindgen(hist[ih])]
      outtruth = otruth_base+sfx2
      outobs = oobs_base+sfx2
      outdeep = odeep_base+sfx2
      outjohnson = ojohnson_base+sfx2
      outflamex = oflamex_base+sfx2
      outndwfs = ondwfs_base+sfx2
      outwise = owise_base+sfx2
      outrcs = orcs_base+sfx2
      outcfhtls = ocfhtls_base+sfx2
      outeuclid = oeuclid_base+sfx2
      outirac = oirac_base+sfx2
      outwise = owise_base+sfx2
      outhsc = ohsc_base+sfx2
      outlsst = olsst_base+sfx2
      outwfirst = owfirst_base+sfx2
      outtwomass = otwomass_base+sfx2
      outvhs = ovhs_base+sfx2
      outviking = oviking_base+sfx2
      outsdss = osdss_base+sfx2
      outstripe82 = ostripe82_base+sfx2
      outbcs = obcs_base+sfx2

      ;;;is this the first file?  If so, we make new output file, if not will append
      if (j eq 0) then begin
        mwrfits, truth[ii], outtruth, /create
        mwrfits, obs[ii], outobs, /create
        mwrfits, deep[ii], outdeep, /create
        mwrfits, johnson[ii], outjohnson, /create
        mwrfits, flamex_out[ii], outflamex, /create
        mwrfits, ndwfs[ii], outndwfs, /create
        mwrfits, wise[ii], outwise, /create
        mwrfits, rcs[ii], outrcs, /create
        mwrfits, cfhtls[ii], outcfhtls, /create
        mwrfits, euclid[ii], outeuclid, /create
        mwrfits, irac[ii], outirac, /create
        mwrfits, wise[ii], outwise, /create
        mwrfits, hsc[ii], outhsc, /create
        mwrfits, lsst[ii], outlsst, /create
        mwrfits, wfirst[ii], outwfirst, /create
        mwrfits, twomass[ii], outtwomass, /create
        mwrfits, vhs[ii], outvhs, /create
        mwrfits, viking[ii], outviking, /create
        mwrfits, sdss[ii], outsdss, /create
	mwrfits, stripe82[ii], outstripe82, /create
	mwrfits, bcs[ii], outbcs, /create
      endif else begin
        ttruth = mrdfits(outtruth,1)
	tobs = mrdfits(outobs,1)
	tdeep = mrdfits(outdeep,1)
	tjohnson = mrdfits(outjohnson,1)
	tflamex = mrdfits(outflamex,1)
	tndwfs = mrdfits(outndwfs,1)
	twise = mrdfits(outwise,1)
	trcs = mrdfits(outrcs,1)
	tcfhtls = mrdfits(outcfhtls,1)
	teuclid = mrdfits(outeuclid,1)
	tirac = mrdfits(outirac,1)
	twise = mrdfits(outwise,1)
	thsc = mrdfits(outhsc,1)
	tlsst = mrdfits(outlsst,1)
	twfirst = mrdfits(outwfirst,1)
	ttwomass = mrdfits(outtwomass,1)
	tvhs = mrdfits(outvhs,1)
	tviking = mrdfits(outviking,1)
	tsdss = mrdfits(outsdss,1)
	tstripe82 = mrdfits(outstripe82, 1)
	tbcs = mrdfits(outbcs, 1)

        mwrfits, [ttruth, truth[ii]], outtruth, /create
        mwrfits, [tobs, obs[ii]], outobs, /create
        mwrfits, [tdeep, deep[ii]], outdeep, /create
        mwrfits, [tjohnson, johnson[ii]], outjohnson, /create
        mwrfits, [tflamex, flamex_out[ii]], outflamex, /create
        mwrfits, [tndwfs, ndwfs[ii]], outndwfs, /create
        mwrfits, [twise, wise[ii]], outwise, /create
        mwrfits, [trcs, rcs[ii]], outrcs, /create
        mwrfits, [tcfhtls, cfhtls[ii]], outcfhtls, /create
        mwrfits, [teuclid, euclid[ii]], outeuclid, /create
        mwrfits, [tirac, irac[ii]], outirac, /create
        mwrfits, [twise, wise[ii]], outwise, /create
        mwrfits, [thsc, hsc[ii]], outhsc, /create
        mwrfits, [tlsst, lsst[ii]], outlsst, /create
        mwrfits, [twfirst, wfirst[ii]], outwfirst, /create
        mwrfits, [ttwomass, twomass[ii]], outtwomass, /create
        mwrfits, [tvhs, vhs[ii]], outvhs, /create
        mwrfits, [tviking, viking[ii]], outviking, /create
        mwrfits, [tsdss, sdss[ii]], outsdss, /create
	mwrfits, [tstripe82, stripe82[ii]], outstripe82, /create
	mwrfits, [tbcs, bcs[ii]], outbcs, /create
        ttruth = 0
        tobs = 0
        tdeep = 0
        tjohnson = 0
        tflamex = 0
        tndwfs = 0
        twise = 0
        trcs = 0
        tcfhtls = 0
        teuclid = 0
        tirac = 0
        twise = 0
        thsc = 0
        tlsst = 0
        twfirst = 0
        ttwomass = 0
        tvhs = 0
        tviking = 0
        tsdss = 0
        tstripe82 = 0
        tbcs = 0

      endelse
    endfor ;;;ends loop around output pixels

    ;;;just delete some data to save memory
    truth = 0
    obs = 0
    deep = 0
    johnson = 0
    flamex = 0
    ndwfs = 0
    wise = 0
    rcs = 0
    cfhtls = 0
    euclid = 0
    irac = 0
    wise = 0
    hsc = 0
    lsst = 0
    wfirst = 0
    twomass = 0
    vhs = 0
    viking = 0
    sdss = 0
    stripe82 = 0
    bcs = 0
  endfor


print, "Finished processing pixel ", pixel

end

