
; specify our input files
dir = '/nfs/slac/g/ki/ki22/cosmo/jderose/addgals/catalogs/Fox/Catalog_v1.0/individual_box_files/'
fbase = 'PO_Fox-1_1050_truth_no_photoz.'
magbase = 'PO_Fox-1_1050_DR8_Mock_mag.'
indices = ['0', '1', '4', '5', '6', '7']

; where we'll save the catalog
outdir = '/nfs/slac/g/ki/ki22/cosmo/jderose/addgals/catalogs/Fox/Catalog_v1.0/wp_cat'
outbase = 'galaxies_wp_cuts.fit'
fout = outdir+outbase

; define our cuts
cspeed = 2.99792e5
mrmin = ['18', '19', '20', '21']
zmax = [12500.0/cspeed, 19250.0/cspeed, 31900.0/cspeed, 47650./cspeed]

; loop over the indices to readin the catalogs, make the temp files
for i = 0, N_ELEMENTS(indices) - 1 do begin

  ; read the data
  g = mrdfits(dir+fbase+indices[i]+'.fit', 1)
  mag = mrdfits(dir+fbase+indices[i]+'.fit', 1)

  ; replace the DES mags with the SDSS ones
  g.amag = mag.amag
  g.omag = mag.omag

  ; cut the file to the depth
  ii = where(g.z le zmax[3] and mag.amag[2] le -1*float(mrmin[0]))

  ; create a new file the first pass
  if (i eq 0) then begin
    mwrfits, g[ii], fout, /create
  ; append to the file on the next pass
  endif else begin
    tg = mrdfits(fout, 1)
    mwrfits, [tg, g[ii]], fout, /create      
  endelse

endfor

end

