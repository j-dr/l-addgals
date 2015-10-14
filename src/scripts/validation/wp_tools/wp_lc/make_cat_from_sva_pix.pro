pro make_cat_from_sva_pix
; specify our input files
dir = '/nfs/slac/g/ki/ki23/des/jderose/addgals/catalogs/chinchilla-tuning/SVA/hpix32/'
fbase = 'Chinchilla-tuning_galaxies.'
magbase = 'Chinchilla-tuning_sdss25.'
indices = ['10074', '10905', '10683', '10448', '11202', '10206', '10907',$
 '10687', '10080', '10456', '11296', '10327', '10903', '10573', '10078',$
 '10452', '11200', '10333', '10801', '10575', '10200', '10454', '11378',$
 '10331', '10799', '10795', '10198', '10569', '11294', '10070', '10689',$
 '10909', '10329', '11012', '10571', '11206', '10072', '10691', '10901',$
 '11105', '10202', '10567', '11204', '10076', '11008', '10797', '10685',$
 '10325', '10450', '11006', '11292', '11107', '10204', '11010', '11111',$
 '11109', '11462', '11464', '11384', '11382', '11466', '11380', '10205',$
 '10207', '10326', '10330', '10332', '10203', '11203', '10334', '10328',$
 '10568', '10453', '11295', '11291', '11205', '10570', '10684', '11293',$
 '11379', '11297', '10451', '10457', '10449', '10576', '10572', '11201',$
 '10455', '10796', '10794', '10686', '10688', '10690', '10574', '10199',$
 '10075', '10197', '10079', '10201', '10073', '10071', '10077', '10798',$
 '10802', '11005', '10904', '10902', '10908', '10906', '10800', '11465',$
 '11461', '11467', '11383', '11381', '11463', '11106', '11009', '11011',$
 '11104', '11110', '11108', '11199', '11007']

; where we'll save the catalog
outdir = '/nfs/slac/g/ki/ki23/des/jderose/addgals/catalogs/chinchilla-tuning/SVA/validation/wp_cat'
outbase = 'galaxies_wp_cuts.fit'
fout = outdir+outbase

; define our cuts
cspeed = 2.99792e5
mrmin = ['18', '19', '20', '21']
zmax = [12500.0/cspeed, 19250.0/cspeed, 31900.0/cspeed, 47650./cspeed]

; loop over the indices to readin the catalogs, make the temp files
for i = 0, N_ELEMENTS(indices) - 1 do begin

  ; read the data
  g = mrdfits(dir+'/'+indices[i]+'/'+fbase+indices[i]+'.fit', 1)
  mag = mrdfits(dir+'/'+indices[i]+'/'+fbase+indices[i]+'.fit', 1)

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

