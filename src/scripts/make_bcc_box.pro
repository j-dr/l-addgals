pro make_bcc_box, dir=dir, $ 
	sim_zmin=sim_zmin, sim_zmax=sim_zmax, $
	nproc=nproc, $
	omegam=omegam, omegal=omegal, $
	simfile=simfile, rnnfile=rnnfile, $
	halofile=halofile, rnn_halofile=rnn_halofile, $
	simname=simname, boxsize=boxsize, denspdfstr=denspdfstr,$
	npix=npix, bcg_mass_lim = bcg_mass_lim, paramfile=paramfile, $
	catdir=catdir, hfile=hfile


;;;parameters that need to be specified
;dir = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Lb1050_v1.0/'
;sim_zmin = 0.0
;sim_zmax = 0.36
;nproc = 8
;omegam = 0.23
;omegal = 0.77
;simname = 'Chinchilla'
;boxsize = '1050'
;simfile = ''
;rnnfile = ''
;halofile = ''
;rnn_halofile = ''
dir += '/'

;;;some intialization steps
spawn, 'mkdir -p '+dir
min_dz = 0.02     ;; smallest delta_z (before buffer) of a processor
buffer_size = 50. ;;overlap region in Mpc
readcol, '~mbusha/projects/development/healpix_cells_in_quartant', pixnum, format = 'l'
if not KEYWORD_SET(npix) then npix = N_ELEMENTS(pixnum)


;;;polynomial to tell us how deep we need to go as a function of z
coefs = [-10.0157, -24.6315, 44.2633, -47.2943, 24.9187, -4.99625]
coefs[0] += 0.2
ncoefs = N_ELEMENTS(coefs)

;;;sample LF for doing depth calculation
blan = [0.0168,-1.03, -20.41]
lmin = -25.0
lmax = -10.
bin = 0.001
nbins = floor((lmax - lmin)/bin)
mag_bin = fltarr(nbins)
for i = 0L, nbins - 1 do begin
  mag_bin(i) = lmin + i*bin
endfor
LF = schechter_mag(mag_bin, blan)

;;;calculate how deep we need to go
nzbins = 2000
zarr = fltarr(nzbins)
ng = fltarr(nzbins)
ng_bin = fltarr(nzbins)
delz = (sim_zmax+0.2-sim_zmin)/nzbins
zarr = sim_zmin + findgen(nzbins)*delz
table = generate_z_of_r_table(omegam, omegal, zmax = sim_zmin+nzbins*delz + 0.2)
rarr = r_of_z(zarr, table)
ng(0) = 0
for i = 1, nzbins - 1 do begin
   this_vol = 4./3*!PI*(rarr(i)^3 - rarr(i-1)^3)
   this_mmin = coefs(0) + coefs(1)*zarr[i] + coefs(2)*zarr[i]^2 + coefs(3)*zarr[i]^3 + coefs(4)*zarr[i]^4 + coefs(5)*zarr[i]^5
   this_mmin = this_mmin + 1.3*(zarr(i) - 0.1) ;;we need to de-evolve this
   if (this_mmin lt -19.6) then this_mmin = -19.6
   if (this_mmin gt -11.0) then this_mmin = -11.0
   ii = where(mag_bin le this_mmin)
   ng(i) = ng(i-1) + this_vol*INT_TABULATED(mag_bin[ii], LF[ii])
   ng_bin(i) = this_vol*INT_TABULATED(mag_bin[ii], LF[ii])
endfor

zmin = fltarr(nproc)
zmax = fltarr(nproc)
mmin = fltarr(nproc)
zmin[0] = sim_zmin
zmax[nproc-1] = sim_zmax
ka = min(abs(zarr - sim_zmax), loc)
ngpp = ng(loc)/nproc
next_ng = ngpp

this_i = nzbins-1
for ip = nproc - 1, 1, -1 do begin
   ;calculate ng/proc to equally divide the rest of the galaxies
   mean_ng = ng(this_i) / (ip+1.) 

   ;how low must we go to not violate our min(delta z) requirement?
   zmin_max = zmax[ip] - min_dz

   ;find where that puts zmin
   this_ng = 0
   while(this_ng lt mean_ng or zarr(this_i) gt zmin_max) do begin
      this_ng += ng_bin(this_i)
      this_i--
   endwhile
   zmin(ip) = zarr(this_i)
   zmax(ip-1) = zarr(this_i)
endfor
;;now we go through and add a buffer zone
for ip = 0, nproc - 1 do begin
   rmin = max([lookbackdist(zmin(ip), omegam, omegal) - buffer_size, 0.])
   rmax = lookbackdist(zmax(ip), omegam, omegal) + buffer_size
   zmin(ip) = max([z_of_r(rmin, table), sim_zmin])
   zmax(ip) = min([z_of_r(rmax, table), sim_zmax])
endfor

LDnum = strarr(nproc)
ldarr = ['60', '65', '68', '70', '72', '75', '78', '80', '82', '85', '88', '90', '92', '95', '97', '98', '99']
;;calcualte the lasdamas snapshot corresponding to this z-range
for i = 0L, nproc - 1 do begin
   tnum = 0
   while (lasdamas_z(ldarr[tnum]) gt zmin[i]) do tnum++
;   LDnum[i] = '000'+strcompress(string(tnum), /remove_all)
;   LDnum[i] = strmid(LDnum[i], strlen(LDnum[i])-3,3)
   LDnum[i] = ldarr[tnum]
endfor

mmin = coefs(0) + coefs(1)*zmin + coefs(2)*zmin^2 + coefs(3)*zmin^3 + coefs(4)*zmin^4 + coefs(5)*zmin^5
mmin = mmin + 1.3*(zmin - 0.1) ;;we need to de-evolve this
ii = where(mmin lt -19.6, count)
if (count gt 0) then mmin[ii] = -19.6
ii = where(mmin gt -11.0, count)
if (count gt 0) then mmin[ii] = -11.0

for i = 0, nproc - 1 do begin
   print, ' '
   print, i
   print, ' '

   ;;define our string variables
   zminstr = ' '+strcompress(string(zmin[i]), /remove_all)
   zmaxstr = ' '+strcompress(string(zmax[i]), /remove_all)
   numstr = ' '+strcompress(string(i),/remove_all)
   if (i lt 100) then numstr = ' 0'+strcompress(string(i), /remove_all)
   if (i lt 10) then numstr = ' 00'+strcompress(string(i), /remove_all)
   num2str = strcompress(string(i),/remove_all)
   if (i lt 100) then num2str = '0'+strcompress(string(i), /remove_all)
   if (i lt 10) then num2str = '00'+strcompress(string(i), /remove_all)

   ;;;define our LF parameters
   mminstr = ' '+strcompress(string(mmin(i)), /remove_all)
   params = get_ages_dsg_params(0.5*(zmin[i]+zmax[i]), /remove_passive_evolution)
   phistar = params[0]
   phistarstr = ' '+strcompress(string(phistar), /remove_all)

   ;;;define out denspdf and lbcg files
   if not KEYWORD_SET(denspdfstr) then denspdfstr = ' /u/ki/mbusha/projects/modules/idl/addgals/rdel/denspdf_Consuelo02_AGES_rescale_0'+LDnum[i]+'.dat'
   denspdfstr = ' '+denspdfstr
   lbcgstr = ' /u/ki/mbusha/projects/addgals/addgalsPO/denspdf/lbcg_relation_Consuelo02_AGES_rescale_0'+LDnum[i]+'.txt'

   ;;;write our LF for this redshift to a temporary file
   tlf_file = 'tLF.dat'
   write_rescaled_lf, tlf_file, params

   ;;;loop over our healpix numbers
   for ipix = 0, npix - 1 do begin
     pixstr = strcompress(string(pixnum(ipix)), /remove_all)

     ;;define our output directories
     this_dir = dir+pixstr+'/'+num2str+'/'
     pathstr = ' '+this_dir+'/hv_output/'
     out_lf_file = this_dir+'/LF.dat'
     print, this_dir
     ;;generate our executables, etc
; Joe removed paramfile for faber tests
;     command = './make_params_files.sh '+zminstr+zmaxstr+pathstr+mminstr+$
;	denspdfstr+lbcgstr+phistarstr+' '+$
;	pixstr+' '+num2str+' '+$
;	this_dir+' '+simfile+' '+' '+rnnfile+' '+$
;	halofile+' '+rnn_halofile+' '+$
;	simname+' '+boxsize+' '+paramfile+' '+bcg_mass_lim
     command = './make_params_files.sh '+zminstr+zmaxstr+pathstr+mminstr+$
	denspdfstr+lbcgstr+phistarstr+' '+$
	pixstr+' '+num2str+' '+$
	this_dir+' '+simfile+' '+' '+rnnfile+' '+$
	halofile+' '+rnn_halofile+' '+$
	simname+' '+boxsize+' '+bcg_mass_lim

     spawn, command
;     spawn, 'cp '+tlf_file+' '+out_lf_file
;     spawn, 'mkdir -p '+this_dir+'/idl'
;     spawn, 'cp hv '+this_dir+'/hv'
;     spawn, 'cp StringParameters '+this_dir
;     spawn, 'cp NumericalParameters '+this_dir

     ;;;seup idl programs
;     idl_file = this_dir+'/idl/run.idl'
;     go_file = this_dir+'/idl/go.sh'
;     spawn, 'echo ".compile ./add_bcgs.pro" > '+idl_file
;     spawn, 'echo "name = '+"'"+'PO_Aardvark_1050_'+pixstr+'.'+num2str+"'"+'"'+'>> '+idl_file
;     spawn, 'echo "create_catalog, name, g, h, /des,/vista,/deep,/johnson,/flamex" >> '+idl_file
;     spawn, 'echo "exit" >> '+idl_file
;     spawn, 'echo "/afs/slac.stanford.edu/g/ek/rsi/idl_6.3/bin/idl run.idl" > '+go_file
;     spawn, 'chmod 744 '+go_file

   endfor
endfor

;;;make the checking/resubmission scripts
cmd = './make_submission_files.sh '+simname+' '+boxsize+' '+dir
spawn, cmd

; make the combining files
fbase = 'PO_'+simname+'_'+boxsize
spawn, 'mkdir -p '+dir+'/make_bcc_file_logs/'
command = './prepare_make_bcc_file.sh '+dir+' '+$
 	fbase+' '+catdir+'/individual_box_files '+hfile+' '+string(nproc)
spawn, command

end
