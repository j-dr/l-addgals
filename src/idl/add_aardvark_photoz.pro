pro add_aardvark_photoz,path, filebase, training_file, weights_file

spawn, 'source /afs/slac/g/ki/root/bin/thisroot.csh'

mag_lim = [24.6, 24.1, 24.4, 23.8, 21.3]

print, "Reading first set of catalogs...."
truth_files = list_with_path(filebase+'_truth_no_photoz.*.fit', path)
visible_files = list_with_path(filebase+'_no_photoz.*.fit', path)
nfiles = N_ELEMENTS(truth_files)/2 ;;not sure why this is detecting files twice

for ifile = 10L, nfiles - 1 do begin
  obs_full = mrdfits(truth_files[ifile], 1, truth_hdr)
  obs_visible = mrdfits(visible_files[ifile], 1, obs_hdr)
  truth_hdr2 = truth_hdr
  sxdelpar, truth_hdr2, 'Z_LO'
  sxdelpar, truth_hdr2, 'Z_HI'
  sxdelpar, truth_hdr2, 'AREA'
  sxaddpar, truth_hdr2, 'AREA', 53.7
  obs_hdr2 = truth_hdr2
  sxaddpar, truth_hdr2, 'Z_LO', 0.0
  sxaddpar, truth_hdr2, 'Z_HI', 2.0


  ;;;extract the healpix number from the filename
  split_pts = strsplit(truth_files[ifile], '.')
  nsplit = N_ELEMENTS(split_pts)
  healpix_num = strmid(truth_files[ifile], split_pts[nsplit-2], $
	               split_pts[nsplit-1]-split_pts[nsplit-2]-1)

  print, "Doing zCarlos and writing sub-regions..."

;     sxdelpar, truth_hdr2, 'RA_MIN'
;     sxdelpar, truth_hdr2, 'RA_MAX'
;     sxaddpar, truth_hdr2, 'RA_MIN', tramin
;     sxaddpar, truth_hdr2, 'RA_MAX', tramax  
;     sxdelpar, obs_hdr2, 'RA_MIN'
;     sxdelpar, obs_hdr2, 'RA_MAX'
;     sxaddpar, obs_hdr2, 'RA_MIN', tramin
;     sxaddpar, obs_hdr2, 'RA_MAX', tramax
  sxaddpar, truth_hdr2, 'HEALPIX_NUM', healpix_num
  sxaddpar, obs_hdr2, 'HEALPIX_NUM', healpix_num

  ;;;run zCarlos
  zcarlos_file =  path+filebase+'_zCarlos.'+healpix_num+'.fit'
  make_carlos_photoz, training_file, truth_files[ifile], $
                      ofile =zcarlos_file, nne = 100
  zCarlos = mrdfits(zcarlos_file,1)
  obs_full.zCarlos = zCarlos.zCarlos_photoz
  obs_visible.zCarlos = zCarlos.zCarlos_photoz

  ;;;run annz
  annz_file =path+filebase+'_annz.'+healpix_num+'.fit'
  run_annz_photoz, truth_files[ifile], annz_file, weights_file
  annz = mrdfits(annz_file,1)
  obs_full.annz = annz.annz_photoz
  obs_full.annz_err = annz.annz_err
  obs_visible.annz = annz.annz_photoz
  obs_visible.annz_err = annz.annz_err

if(0) then begin
  ;;;run arborz
  arborz_file = path+filebase+'_arborz.'+healpix_num+'.fit'
  run_arborz_photoz, truth_files[ifile], outfile=arborz_file
  arborz = mrdfits(arborz_file,1)
  obs_full.arborz = arborz.arborz_photoz
  obs_full.arborz_err = arborz.arborz_error
  obs_visible.arborz = arborz.arborz_photoz
  obs_visible.arborz_err = arborz.arborz_error
endif

  ;;;write out our sectional files
  mwrfits, obs_visible, path+filebase+'.'+healpix_num+'.fit', $
	   obs_hdr2, /create
  mwrfits, obs_full, path+filebase+'_truth.'+healpix_num+'.fit', $
	   truth_hdr2, /create
endfor

end



