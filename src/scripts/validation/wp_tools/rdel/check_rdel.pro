;pro check_rdel, path, denspdf_file, mag=mag, sim_file = sim_file, $
pro check_rdel, mag=mag, carr=carr, $
		mr=mr, rdel=rdel, zred=zred, $
                hv_path=hv_path, $
                gfile=gfile, $
                shamfile=shamfile, sham_mr=sham_mr, sham_rdel=sham_rdel, $
		paramfile=paramfile, $
                zrange=zrange, simtype=simtype, Q=Q

if not KEYWORD_SET(simtype) then begin
  print, "Error!  Simtype not specified.  "
  return
endif

; setup the samples to be plotted
if not KEYWORD_SET(mag) then mag = [-18, -19, -20, -21, -22]
if not KEYWORD_SET(carr) then carr = [!purple, !blue, !green, !yellow, !red]
nmag = N_ELEMENTS(mag)

; read in the data -- could be different formats
enough_inputs = 0
print, "Reading addgals data..."
if (keyword_set(mr) and keyword_set(rdel)) then begin
  if (simtype eq 'LC' and not KEYWORD_SET(zred)) then begin
    print, "Error!  Want to do a LC simulation but didn't specify zred!"
    return
  endif
  if (N_ELEMENTS(mr) eq N_ELEMENTS(rdel)) then begin
    enough_inputs = 1
  endif else begin
    print, "Error!  Inconsistent mr and rdel specified!"
    print, "N_ELEMENTS(mg) = ", N_ELEMENTS(mr)
    print, "N_ELEMENTS(rdel) = ", N_ELEMENTS(rdel)
    return
  endelse
endif else if (KEYWORD_SET(hv_path)) then begin
  readcol, hv_path+'/gal_ginfo1.dat', sedid, mr, ra, dec, zred, cent, id, $
	format='l,f,f,f,f,l,l'
  rdfloat, hv_path+'/gal_dinfo.dat', grnn, rnn
  ;rdel = rnn
  rdel = grnn
  enough_inputs = 1
endif else if (KEYWORD_SET(gfile)) then begin
  g = mrdfits(gfile, 1)
  mr = g.amag[2]
  rdel = g.rnn
  zred = g.z
  enough_inputs = 1
endif else begin
  print, "Error!  Didn't not specify enough inputs to read the modeled data!"
  return
endelse

; remove passive evolution
if KEYWORD_SET(Q) then begin
  print, "Removing passive evolution."
  mr -= Q*(1./(1+zred) - 1./1.1)
endif

; read the model fit, if requested
if KEYWORD_SET(paramfile) then begin
  print, "reading parameters..."
  cmparams = get_cm_params(paramfile)
  csparams = get_cs_params(paramfile)
  fmparams = get_fm_params(paramfile)
  fsparams = get_fs_params(paramfile)
  pparams = get_p_params(paramfile)
endif

; read simdata, if requested
if KEYWORD_SET(shamfile) then begin
  print, "reading SHAM data..."
  if not (KEYWORD_SET(sham_mr) and KEYWORD_SET(sham_rdel)) then $
    readcol, shamfile, sham_id, sham_mr, sham_rdel, format='l,f,f'
endif

; setup our plotting array
rarr = (1+findgen(1000))*0.01

; cut our sample data
if KEYWORD_SET(zrange) then begin
   ii = where(zred ge zrange[0] and zred le zrange[1])
   mr = mr[ii]
   rdel=rdel[ii]
   zred = zred[ii]
endif

; we're ready to start plotting!
print, "making plots..."
for i = 0, nmag - 1 do begin

  ; setup an empty plot
  title = 'Mr < '+strcompress(string(mag[i]), /remove_all)
  plot, [0,0], [0,0], xrange = [0,8], yrange = [0,1], /nodat, $
	xtitle = 'r!Ddel', title=title

  ; cut our sample and plot it
  ii = where(mr le mag[i])
  if (N_ELEMENTS(ii) lt 2) then continue
  pplothist, rdel[ii], bin = 0.05, /prob, /over

  ; add the model
  if KEYWORD_SET(paramfile) then begin
    zmid = mean(zred[ii])
    cm = cm_fcn3(mag[i], zmid, cmparams)
    cs = cs_fcn3(mag[i], zmid, csparams)
    fm = fm_fcn3(mag[i], zmid, fmparams)
    fs = fs_fcn3(mag[i], zmid, fsparams)
    p = p_fcn3(mag[i], zmid, pparams)

    oplot, rarr, rdel_fcn(rarr, [cm, cs, fm, fs, p]), linestyle = 1, color=!red
  endif

  ; add the sham data
  if KEYWORD_SET(shamfile) then begin
    ii = where(sham_mr le mag[i])
    pplothist, sham_rdel(ii), bin = 0.03, /prob, /over, linestyle = 2, color=!blue
  endif

  if (mag[i] eq -19) then begin
    oplot, rarr, rdel_fcn(rarr, [-0.467791, 0.510618, 2.88317, 1.12333, 0.726224]), color = !green
    print, cm, cs, fm, fs, p
  endif

  ; the all-important legend
  legend, ['ADDGALS', 'Model', 'SHAM'], linestyle = [0,1,2], /right
endfor

stop

end
