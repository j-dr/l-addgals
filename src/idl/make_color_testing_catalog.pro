;;;dir = '~/ki01/projects/addgals/Carmen/newcode/v1/'

nproc = 10
;;;dir = '/nfs/slac/g/ki/ki11/des/mbusha/catalogs/Carmen/LC/v4.40/v13/'
;zmin = [0.0, 0.23, 0.35, 0.45, 0.55, 0.64, 0.74, 0.85, 0.98, 1.13]
;zmax = [0.23, 0.35, 0.45, 0.55, 0.64, 0.74, 0.85, 0.98, 1.13, 1.33]

dir = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchilla-tuning/analysis/lc/dc5/addgals/'
zmin = [0.0, 0.69, 0.91, 1.09, 1.25, 1.39, 1.53, 1.66, 1.79, 1.91]
zmax = [0.69, 0.91, 1.09, 1.25, 1.39, 1.53, 1.66, 1.79, 1.91, 2.0]

spawn, 'mkdir '+dir



;;;loop through processors and copy the relevant files
for i = 0, nproc - 1 do begin
   ;;;output diagnostics to see what slice we're on
   print, ' '
   print, i
   print, ' '

   ;;;change some variables to strings
   zminstr = ' '+strcompress(string(zmin[i]), /remove_all)
   zmaxstr = ' '+strcompress(string(zmax[i]), /remove_all)
   numstr = ' '+strcompress(string(i),/remove_all)
   if (i lt 100) then numstr = ' 0'+strcompress(string(i), /remove_all)
   if (i lt 10) then numstr = ' 00'+strcompress(string(i), /remove_all)
   num2str = strcompress(string(i),/remove_all)
   if (i lt 100) then num2str = '0'+strcompress(string(i), /remove_all)
   if (i lt 10) then num2str = '00'+strcompress(string(i), /remove_all)
   mminstr = ' -11.0' ;;isn't used for color assignment
   pathstr = ' '+dir+num2str+'/hv_output/'
   spawn, 'mkdir '+dir+num2str

   ;;;denspdf and lbcg files aren't used, but they still need to be defined
   denspdfstr = ' /u/ki/mbusha/projects/modules/idl/addgals/rdel/denspdf_Consuelo_scatter0.20_histograms_individual_params_099.txt'
   lbcgstr = ' /afs/slac.stanford.edu/u/ki/mbusha/projects/modules/idl/addgals/rdel/lbcg_m200_Consuelo_scatter0.20_histograms_099.txt'

   ;phistar needs to be defined, even if it isn't used 
   params = get_ages_dsg_params(0.5*(zmin[i]+zmax[i]), /remove_passive_evolution)
   phistar = params[0]
   phistarstr = ' '+strcompress(string(phistar), /remove_all)

   ;;;create parameter files based on everytyhing we've wpecified (really just redshift ranges)
   command = './make_params_files_220sqdeg.sh '+zminstr+zmaxstr+pathstr+mminstr+denspdfstr+lbcgstr+phistarstr
   spawn, command
   spawn, 'mkdir '+dir+num2str+'/idl'
   spawn, 'cp hv '+dir+num2str+'/hv.'+num2str
   spawn, 'cp StringParameters '+dir+num2str
   spawn, 'cp NumericalParameters '+dir+num2str
   idl_file = dir+num2str+'/idl/run.idl'
   go_file = dir+num2str+'/idl/go.sh'
   spawn, 'echo "cp /afs/slac.stanford.edu/u/ki/mbusha/projects/addgals/idl/make_catalog.sav ." > '+go_file
   spawn, 'echo "/afs/slac.stanford.edu/u/ki/mbusha/bin/idl_vm_run.py make_catalog.sav PO_Carmen02_000" >> '+go_file
   spawn, 'chmod 744 '+go_file

   ;;;link the galaxy .ascii file to the relevant directory
   ;galaxies_file = '/nfs/slac/g/ki/ki11/des/mbusha/catalogs/Carmen/LC/v4.40/v08/DES_catalog/galaxies'+strcompress(string(i),/remove_all)+'.ascii'
   galaxies_file = '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchilla-tuning/analysis/lc/dc5/galaxies'+strcompress(string(i),/remove_all)+'.ascii'
   spawn, 'ln -s '+galaxies_file+' '+dir+num2str+'/galaxies.ascii'

endfor

spawn, 'cp runscript ' + dir
spawn, 'cp run.sh ' + dir
spawn, 'cp submit_jobs ' + dir

end
