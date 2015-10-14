pro calculate_full_wp, box_size, outpath, $
	hv_path=hv_path, $
	rockfile=rockfile, shamfile=shamfile

spawn, 'mkdir -p '+outpath

; read the input data
if KEYWORD_SET(hv_path) then begin
  rdfloat, hv_path+'/gal_pinfo.dat', px, py, pz, vx, vy, vz
  readcol, hv_path+'/gal_ginfo1.dat', sedid, mr, ra, dec, zred, cent, id, $
        format='l,f,f,f,f,l,l'
endif else if (KEYWORD_SET(rockfile) and KEYWORD_SET(shamfile)) then begin
  read_hlist_general, rockfile, h
  px = h.x
  py = h.y
  pz = h.z
  vx = h.vx
  vy = h.vy
  vz = h.vz
  readcol, shamfile, id, mr, rdel, format = 'll,f,f'
endif else begin
  print, "Error!  Not enough inputs specified!"
  print, "You must specify either hv_path OR rockfile and shamfile"
  return
endelse

; set the parameters for the limits of the calculation
mag_thresh = [-18, -19, -20, -21]
PIMAX = ['40', '40', '80', '80']
nsamples = N_ELEMENTS(mag_thresh)

; executable where we'll store commands for generating the wp calculation
makewp_file = outpath+'makewp_cameron.sh'
spawn, 'rm -f '+makewp_file
spawn, 'touch '+makewp_file
spawn, 'chmod 744 '+makewp_file

; loop over all the samples we're trying to calcualte wp for
for i = 0, nsamples - 1 do begin

  ; setup the file for storing the wp data
  num = strcompress(string(abs(mag_thresh[i])), /remove_all)
  outfile = outpath+'/gal_wpinfo_mt'+num+'_all.dat'

  ; select the galaxies we'll calculate on
  ii = where(mr le mag_thresh[i]) 

  ; make a reduced file.  
  openw,1,outfile
  niceprintf,1,px[ii], py[ii], pz[ii], vx[ii], vy[ii], vz[ii]
  close, 1

  ; run the calculation 
  cmd = 'echo "~mbusha/projects/jeremy_wp/wp_covar 0.1 20 15 '+string(box_size)+' 0 '+string(box_size)+' 1 '+outfile+' a 1 1 1  > '+outpath+'/wp_mt'+num+'_all.dat" >> '+makewp_file
  print, cmd
  spawn, cmd

endfor



end


