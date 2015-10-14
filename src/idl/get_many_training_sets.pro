pro get_many_training_sets, path, fbase, sbase, obase, sq_deg, frac

npatches = 900
nside = 8

ramin = 0.0
ramax = 180.0
decmin = 0.0
decmax = 90.0
thetamax = (90.-decmin)*!PI/180.
thetamin = (90.-decmax)*!PI/180.
costhetamax = cos(thetamax)
costhetamin = cos(thetamin)

n_areas = N_ELEMENTS(obase)

file_old = 'junk'
for i = 0L, npatches - 1 do begin
  ;;readin the file that contains this pixel, if necessary
  strip = strcompress(string(i), /remove_all)
  file = path+'/truth/'+fbase+'.'+strip+'.fit'
  sfile = path+'/DR8/'+sbase+'.'+strip+'.fit'
  if (file_test(file) eq 0) then continue
  print, "doing patch "+strcompress(string(i),/remove_all)
  g = mrdfits(file, 1)
  s = mrdfits(sfile, 1)
  sp = cut_to_dr8_training_limits(s, ii=ii)
  gp = g[ii]

  for ia = 0L, n_areas - 1 do begin

    ;;define a window around the pixel
    dec = mean(g.dec)
    ra = mean(g.ra)
    ddec = sqrt(sq_deg[ia])
    dec1 = dec - 0.5*ddec
    dec2 = dec1 + ddec
    theta1 = (90.-dec1)*!PI/180.
    theta2 = (90.-dec2)*!PI/180.
    dra = (sq_deg[ia]*(!PI/180.)^2/(cos(theta2) - cos(theta1)))*180./!PI
    ra1 = ra - 0.5*dra
    ra2 = ra + 0.5*dra

    ;;identify all galaxies in the window
    ii = where(gp.ra ge ra1 and gp.ra le ra2 and $
  	       gp.dec ge dec1 and gp.dec le dec2, count)
    jj = get_unique_sample(ii, frac[ia])
    print, 'Number of galaxies in patch '+strcompress(string(i),/remove_all)+', model '+strcompress(string(ia),/remove_all)+' : ', count

    ;;;save these galaxies
    ofile = obase[ia]+'_Truth.'+strcompress(string(i),/remove_all)+'.fit'
    sofile = obase[ia]+'_sdss_mag.'+strcompress(string(i),/remove_all)+'.fit'
    mwrfits, gp[ii[jj]], ofile, /create
    mwrfits, sp[ii[jj]], sofile, /create
  endfor

endfor


end

