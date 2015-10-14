pro add_bcc_arborz,path

  spawn, 'source /afs/slac/g/ki/root/bin/thisroot.csh'

  print, "Reading first set of catalogs...."
  obs_full = mrdfits(path+'DES_Mock_Baseline_truth.fit', 1, truth_hdr)
  obs_visible = mrdfits(path+'DES_Mock_Baseline.fit', 1, obs_hdr)
  training_file = path+'DES_Mock_training_set.fit'
  validation_file = path+'DES_Mock_validation_set.fit'
  truth_hdr2 = truth_hdr
  obs_hdr2 = obs_hdr

  print, "Training ArborZ..."
  train_arborz_photoz, training_file

  print, "Doing ArborZ in sub-regions..."
  ramin = 10
  dra = 2
  for i = 0, 9 do begin
     tramin = ramin + dra*i
     tramax = tramin + dra
     ii = where(obs_full.ra ge tramin and obs_full.ra lt tramax)
     help, ii

     gfile = path+'/DES_Mock_Baseline_truth.'+strcompress(string(i),/remove_all)+'.fit'
     sxdelpar, truth_hdr2, 'RA_MIN'
     sxdelpar, truth_hdr2, 'RA_MAX'
     sxaddpar, truth_hdr2, 'RA_MIN', tramin
     sxaddpar, truth_hdr2, 'RA_MAX', tramax  
     sxdelpar, obs_hdr2, 'RA_MIN'
     sxdelpar, obs_hdr2, 'RA_MAX'
     sxaddpar, obs_hdr2, 'RA_MIN', tramin
     sxaddpar, obs_hdr2, 'RA_MAX', tramax

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



