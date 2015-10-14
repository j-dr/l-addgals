;;;dir = '~/ki01/projects/addgals/Carmen/newcode/v7/'
dir = '/nfs/slac/g/ki/ki11/des/mbusha/catalogs/Carmen/Octant/v4.20/v02/'
min_dz = 0.02     ;; smallest delta_z (before buffer) of a processor
buffer_size = 50. ;;overlap region in Mpc

sim_zmin = 0.0
sim_zmax = 0.33
nproc = 16

omegam = 0.25
omegal = 0.75


;coefs = [-9.92999, -46.9112, 142.519, -225.178, 162.726, -43.0895]
;coefs[0] -= 2.5
coefs = [-7.61893, -77.8381, 216.848, -314.270, 222.880, -61.2050] ;calibrated for DR8 depth


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


table = generate_z_of_r_table(omegam, omegal, zmax = sim_zmax + 0.2)
nzbins = 2000
zarr = fltarr(nzbins)
ng = fltarr(nzbins)
ng_bin = fltarr(nzbins)
delz = (sim_zmax+0.2-sim_zmin)/nzbins
zarr = sim_zmin + findgen(nzbins)*delz
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
ldarr = ['060', '065', '068', '070', '072', '075', '078', '080', '082', '085', '088', '090', '092', '095', '097', '098', '099']
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

spawn, 'mkdir -p '+dir
for i = 0, nproc - 1 do begin
   print, ' '
   print, i
   print, ' '
   zminstr = ' '+strcompress(string(zmin[i]), /remove_all)
   zmaxstr = ' '+strcompress(string(zmax[i]), /remove_all)
   numstr = ' '+strcompress(string(i),/remove_all)
   if (i lt 100) then numstr = ' 0'+strcompress(string(i), /remove_all)
   if (i lt 10) then numstr = ' 00'+strcompress(string(i), /remove_all)
   num2str = strcompress(string(i),/remove_all)
   if (i lt 100) then num2str = '0'+strcompress(string(i), /remove_all)
   if (i lt 10) then num2str = '00'+strcompress(string(i), /remove_all)
   spawn, 'mkdir '+dir+num2str
   mminstr = ' '+strcompress(string(mmin(i)), /remove_all)
   pathstr = ' '+dir+num2str+'/hv_output/'
   ;denspdfstr = ' /u/ki/mbusha/projects/addgals/addgalsPO/denspdf/denspdf_Consuelo02_'+LDnum[i]+'_vmax.dat'
   ;lbcgstr = ' /u/ki/mbusha/projects/addgals/addgalsPO/denspdf/lbcg_relation_Consuelo02_'+LDnum[i]+'.txt'
   if (1) then begin ;doing the AGES parameterization
     denspdfstr = ' /u/ki/mbusha/projects/modules/idl/addgals/rdel/denspdf_Consuelo02_AGES_rescale_'+LDnum[i]+'.dat'
;     lbcgstr = ' /u/ki/mbusha/projects/addgals/addgalsPO/denspdf/lbcg_relation_Consuelo02_vir_AGES_'+LDnum[i]+'.txt'
     lbcgstr = ' /u/ki/mbusha/projects/addgals/addgalsPO/denspdf/lbcg_relation_Consuelo02_AGES_rescale_'+LDnum[i]+'.txt'
     params = get_ages_dsg_params(0.5*(zmin[i]+zmax[i]), /remove_passive_evolution)
     ;params = get_ages_params(0.5*(zmin[i]+zmax[i]), /remove_passive_evolution)
     ;params = get_ages_params(0.05, /remove_passive_evolution)
     phistar = params[0]
     out_lf_file = dir+num2str+'/LF.dat'
     write_rescaled_lf, out_lf_file, params
     denspdfstr = dir+num2str+'/denspdf.dat'
     write_start_denspdf_file, denspdfstr, 0.5*(zmin[i]+zmax[i])
   endif
   phistarstr = ' '+strcompress(string(phistar), /remove_all)

   command = './make_params_files.sh '+zminstr+' '+zmaxstr+' '+pathstr+' '+mminstr+' '+denspdfstr+' '+lbcgstr+' '+phistarstr
   spawn, command
   spawn, 'mkdir '+dir+num2str+'/idl'
   spawn, 'cp hv '+dir+num2str+'/hv.'+num2str
   spawn, 'cp StringParameters '+dir+num2str
   spawn, 'cp NumericalParameters '+dir+num2str
   idl_file = dir+num2str+'/idl/run.idl'
   go_file = dir+num2str+'/idl/go.sh'
;   spawn, 'echo ".compile ./add_bcgs.pro" > '+idl_file
;   spawn, 'echo "name = '+"'"+'PO_Carmen02_'+num2str+"'"+'"'+'>> '+idl_file
;   spawn, 'echo "create_catalog, name, g, h, /des,/vista,/deep,/johnson,/flamex" >> '+idl_file
;   spawn, 'echo "exit" >> '+idl_file
;   spawn, 'echo "/afs/slac.stanford.edu/g/ek/rsi/idl_6.3/bin/idl run.idl" > '+go_file
   spawn, 'echo "cp /afs/slac.stanford.edu/u/ki/mbusha/projects/addgals/idl/make_catalog.sav ." > '+go_file
   spawn, 'echo "/afs/slac.stanford.edu/u/ki/mbusha/bin/idl_vm_run.py make_catalog.sav PO_Carmen02_000" >> '+go_file
   spawn, 'chmod 744 '+go_file

endfor

end
