;pro mask_bcc_pixel,path, inbase, outbase, istr
pro mask_bcc_pixel

args=command_line_args(count=ct)
inbase = args[0]
outbase = args[1]
istr = args[2]

;path = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.3/'
;fbase = 'PO_Aardvark'
;istr = '343'


;;;define our file locations
;if not KEYWORD_SET(fbase) then fbase = 'PO_Aardvark'
;inbase = path+'/obs_rotated/'+fbase+'_no_photoz_des_masked.'
;outbase = path+'/mask/'+fbase+'_mask.'

;;;define out mask locations
sbands = ['g', 'r', 'i', 'z', 'Y']
maskfile = '/nfs/slac/g/ki/ki11/des/mswanson/DESdata/BCC/simplebccmask/simplebccmask_round82_t_'+sbands+'.pol'

;;;structure to hold the mask information
m1 = create_struct('ID', 0LL, 'polygon', lonarr(5))

;for i = 0L, 800 - 1 do begin
;  istr = strcompress(string(i),/remove_all)
  infile = inbase+'.'+istr+'.fit'
  outfile = outbase+'.'+istr+'.fit'
  tout = outbase+'_'+sbands+'.'+istr+'.txt'
;  if (file_test(infile) eq 0) then continue
  print, " "
  print, "Masking file: "+infile
  print, " "
  g = mrdfits(infile, 1)
  m = replicate(m1, N_ELEMENTS(g))
  m.id = g.id
  for j = 0, 4 do begin
    mask_galaxies, g, tout[j], mask = maskfile[j], sfx = istr
    readcol, tout[j], tra, tdec, tpoly, format = 'f,f,l'
    m.polygon[j] = tpoly
  endfor
  mwrfits, m, outfile, /create
  for j = 0, 4 do spawn, 'rm -f '+tout[j]
;endfor


end
