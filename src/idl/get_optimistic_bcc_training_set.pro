pro get_optimistic_bcc_training_set, path, fbase, outfile=outfile

if not KEYWORD_SET(outfile) then outfile = path+'/BCC_Optimistic_training_set.fit'

npatches = 150
ngpp = 400
mlim = 24.0
sq_deg = 0.8
nside = 8

ramin = 0.0
ramax = 180.0
decmin = 0.0
decmax = 90.0
thetamax = (90.-decmin)*!PI/180.
thetamin = (90.-decmax)*!PI/180.
costhetamax = cos(thetamax)
costhetamin = cos(thetamin)

;;;select some random locations for our patches
ra = randomu(seed,npatches)*(ramax-ramin) + ramin
phi = ra*!PI/180.
theta = acos(randomu(seed,npatches)*(costhetamax-costhetamin) + costhetamin)
dec = 90-theta*180./!PI
ang2pix_ring, nside, theta, phi, ip

;;;sort them in healpix order to sepeed up the i/o
ipsrt = sort(ip)
theta = theta[ipsrt]
phi = phi[ipsrt]
ip = ip[ipsrt]
ra = ra[ipsrt]
dec = dec[ipsrt]

file_old = 'junk'
for i = 0L, npatches - 1 do begin
  print, "doing patch "+strcompress(string(i),/remove_all)
  ;;readin the file that contains this pixel, if necessary
  strip = strcompress(string(ip[i]), /remove_all)
  file = path+'/'+fbase+'.'+strip+'.fit'
  if (file ne file_old) then g = mrdfits(file, 1)
  file_old = file

  ;;define a window around the pixel
  ddec = sqrt(sq_deg)
  dec1 = dec[i] - 0.5*ddec
  dec2 = dec1 + ddec
  theta1 = (90.-dec1)*!PI/180.
  theta2 = (90.-dec2)*!PI/180.
  dra = (sq_deg*(!PI/180.)^2/(cos(theta2) - cos(theta1)))*180./!PI
  ra1 = ra[i] - 0.5*dra
  ra2 = ra[i] + 0.5*dra

  ;;identify all galaxies in the window
  ii = where(g.ra ge ra1 and g.ra le ra2 and $
	     g.dec ge dec1 and g.dec le dec2 and $
	     g.omag[2] lt mlim, count)
  print, 'Number of galaxies in patch '+strcompress(string(i),/remove_all)+': ', count

  ;;select a unique random subset of these galaxies
  if (count le ngpp) then ind = ii
  ind = ii[randomu(seed,ngpp)*count]
  srt = sort(ind)
  u = uniq(ind, srt)
  while (N_ELEMENTS(u) ne ngpp) do begin
    ind = ii[randomu(seed,ngpp)*count]
    srt = sort(ind)
    u = uniq(ind, srt)
  endwhile
  
  ;;;add the galaxies to our training set
  if (i eq 0) then begin
    gspec = g[ind]
  endif else begin
    gspec = [gspec, g[ind]]
  endelse

endfor

;mwrfits, gspec, path+'/BCC_Optimistic_training_set.fit', /create
mwrfits, gspec, outfile, /create

end

