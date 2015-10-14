pro combine_reformatted_addstars_files, inpath, outpath, outbase

for i = 0L, 800 do begin
  istr = strcompress(string(i),/remove_all)
  tcmd = 'ls '+inpath+'/*_truth.'+istr+'.fit'
  ocmd = 'ls '+inpath+'/*_obs.'+istr+'.fit'
  spawn, tcmd, tfiles
  spawn, tcmd, ofiles
  if (tfiles[0] eq '') then continue
  nfiles_pixel = N_ELEMENTS(tfiles)
  nfiles_pixel_str = strcompress(string(nfiles_pixel), /remove_all)
  print, "reprocessing pixel "+istr+" ("+nfiles_pixel_str+" files)"
  truth = mrdfits(tfiles[0],1)
  obs = mrdfits(ofiles[0],1)
  for j = 1L, nfiles_pixel-1 do begin
    ttruth = mrdfits(tfiles[j],1)
    tobs = mrdfits(ofiles[j],1)
    truth = [truth,ttruth]
    obs = [obs, tobs]
  endfor

  truth_file = outpath+'/'+outbase+'_truth_stars.'+istr+'.fit'
  obs_file = outpath+'/'+outbase+'_stars.'+istr+'.fit'
  
  mwrfits, truth, truth_file, /create
  mwrfits, obs, obs_file, /create
endfor

end
