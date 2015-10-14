pro read_add_des_errors, path, outpath, fbase, healpix_num

  g = mrdfits( path+fbase+healpix_num+'.fit', 1 )
  ng = n_elements( g )
  sdss_mag1 = create_struct('amag', fltarr(5), 'omag', fltarr(5), 'omagerr', fltarr(5), 'flux', fltarr(5), 'ivar', fltarr(5))
  sdss_mag = replicate(sdss_mag1, ng)
  sdss_mag.amag = g.amag
  mock_error_apply, 'DES', sdss_mag.amag, flux, ivar, omag, omagerr
  
  sdss_mag.omag = omag
  sdss_mag.omagerr = omagerr
  sdss_mag.flux = flux
  sdss_mag.ivar = 1./ivar^2
  mwrfits, sdss_mag, outpath+fbase+'_DES_Mock_mag.'+healpix_num+'.fit', /create

end

