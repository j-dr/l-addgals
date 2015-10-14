pro add_bcc_zcarlos,path

  spawn, 'source /afs/slac/g/ki/root/bin/thisroot.csh'

;  mag_lim = [24.6, 24.1, 24.4, 23.8, 21.3]

  print, "Reading first set of catalogs...."
  obs_full = mrdfits(path+'DES_Mock_Baseline_truth_no_photoz.fit', 1, truth_hdr)
  obs_visible = mrdfits(path+'DES_Mock_Baseline_no_photoz.fit', 1, obs_hdr)
  training_file = path+'DES_Mock_training_set.fit'
  validation_file = path+'DES_Mock_validation_set.fit'
  weights_file = path+'DES_Mock_annz.wts'
  truth_hdr2 = truth_hdr
  obs_hdr2 = obs_hdr

;  for i = 0L, 4 do begin
;     good = where(obs_full.omag[i] le mag_lim[i], comp=bad)
;     if (bad[0] ge 0) then obs_full[bad].omag[i] = 99
;  endfor

  obs_visible.mag_g = obs_full.omag[0]
  obs_visible.mag_r = obs_full.omag[1]
  obs_visible.mag_i = obs_full.omag[2]
  obs_visible.mag_z = obs_full.omag[3]
  obs_visible.mag_y = obs_full.omag[4]


;  print, "Getting photo-z training set..."
;  ii = get_training_set(obs_full, random=0, nrandom = 3000, /stripe82, deep=deep, /vvds, /zcosmos)
;  srt = sort(obs_full[ii].ra)
;  validation = ii[srt[0:19999]]
;  training = ii[srt[20000:N_ELEMENTS(ii)-1]]
;  mwrfits, obs_full[training], training_file, /create
;  mwrfits, obs_full[validation], validation_file, /create


  print, "Training ArborZ..."
  train_arborz_photoz, training_file

  print, "Doing zCarlos and writing sub-regions..."
  ;;make sub regions that are just 2 degrees wide
  ramin = 10
  dra = 2
  for i = 0, 9 do begin
     tramin = ramin + dra*i
     tramax = tramin + dra
     print, " "
     print, "Doing region with ra: ", tramin, tramax
     print, " "

     sxdelpar, truth_hdr2, 'RA_MIN'
     sxdelpar, truth_hdr2, 'RA_MAX'
     sxaddpar, truth_hdr2, 'RA_MIN', tramin
     sxaddpar, truth_hdr2, 'RA_MAX', tramax  
     sxdelpar, obs_hdr2, 'RA_MIN'
     sxdelpar, obs_hdr2, 'RA_MAX'
     sxaddpar, obs_hdr2, 'RA_MIN', tramin
     sxaddpar, obs_hdr2, 'RA_MAX', tramax

     ii = where(obs_full.ra ge tramin and obs_full.ra lt tramax)
     gfile = path+'/DES_Mock_Baseline_truth.'+strcompress(string(i),/remove_all)+'.fit'
     annz_file =path+'DES_Mock_annz.'+strcompress(string(i),/remove_all)+'.fit'
     mwrfits, obs_full[ii], gfile, /create

     ;;;run zCarlos
     zcarlos_file =  path+'/DES_Mock_zCarlos.'+strcompress(string(i), /remove_all)+'.fit'
     make_carlos_photoz, path+'/DES_Mock_training_set.fit', gfile, $
                         ofile =zcarlos_file, nne = 100
     zCarlos = mrdfits(zcarlos_file,1)
     obs_full[ii].zCarlos = zCarlos.zCarlos_photoz
     obs_visible[ii].zCarlos = zCarlos.zCarlos_photoz

     ;;;run annz
     run_annz_photoz, gfile, annz_file, weights_file
     annz = mrdfits(annz_file,1)
     obs_full[ii].annz = annz.annz_photoz
     obs_full[ii].annz_err = annz.annz_err
     obs_visible[ii].annz = annz.annz_photoz
     obs_visible[ii].annz_err = annz.annz_err

     ;;;run arborz
     arborz_file = path+'DES_Mock_arborz.'+strcompress(string(i),/remove_all)+'.fit'
     run_arborz_photoz, gfile, outfile=arborz_file
     arborz = mrdfits(arborz_file,1)
     obs_full[ii].arborz = arborz.arborz_photoz
     obs_full[ii].arborz_err = arborz.arborz_error
     obs_visible[ii].arborz = arborz.arborz_photoz
     obs_visible[ii].arborz_err = arborz.arborz_error

     ;;;write out our sectional files
     mwrfits, obs_visible[ii], path+'/DES_Mock_Baseline.'+strcompress(string(i),/remove_all)+'.fit', obs_hdr2, /create
     mwrfits, obs_full[ii], path+'/DES_Mock_Baseline_truth.'+strcompress(string(i),/remove_all)+'.fit', truth_hdr2, /create
  endfor

  mwrfits, obs_full, path+'DES_Mock_Baseline_truth.fit', truth_hdr, /create
  mwrfits, obs_visible, path+'DES_Mock_Baseline.fit', obs_hdr, /create



end



