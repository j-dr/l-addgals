pro calculate_full_wp_sham, rockfile, shamfile, box_size, outpath

spawn, 'mkdir -p '+outpath

; read the input data
read_hlist_general, rockfile, h
readcol, shamfile, id, mr, rdel, format = 'll,f,f'

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
  niceprintf,1,h[ii].x, h[ii].y, h[ii].z, h[ii].vx, h[ii].vy, h[ii].vz
  close, 1

  ; run the calculation 
  cmd = 'echo "~/projects/jeremy_wp/wp_covar 0.1 20 15 '+string(box_size)+' 0 '+string(box_size)+' 1 '+outfile+' a 1 1 1  > '+outpath+'/wp_mt'+num+'_all.dat" >> '+makewp_file
  print, cmd
  spawn, cmd

endfor



end


