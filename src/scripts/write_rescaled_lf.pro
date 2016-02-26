pro write_rescaled_lf, outfile, params

;outfile = 'kaka'
;params = get_ages_dsg_params(0.05, /remove_passive_evolution)
;params = get_ages_dsg_params(0.05)

if (0) then begin
rdfloat, '/afs/slac.stanford.edu/u/ki/mbusha/projects/modules/idl/cosmology/sham/blanton_jt_hybrid.dat', tmag, tphi1, tphi
ii = where(tmag lt -19 and tmag gt -22.5)
start_params = [0.0159, 1.05, -20.58]
new_params = mpfitfun('log_schechter_mag', tmag[ii], alog10(tphi[ii]), 1, start_params, /quiet)
phi_rescale = (params[0]/new_params[0])
mstar_rescale = params[2] - new_params[2]
tmag2 = tmag + mstar_rescale
endif

;;;;;;params = get_double_schechter_gaussian_params()
;tmag2 = lindgen(1500)/100. - 25
tmag2 = lindgen(2500)/100. - 25
tmag = tmag2
tphi = double_schechter_gaussian(tmag2, params)
phi_rescale = 1

openw,1,outfile
niceprintf, 1, tmag2, tphi*phi_rescale
close,1

blan = [1.49e-2,-1.05, -20.44]
old_params = get_double_schechter_gaussian_params()

;plot, tmag, double_schechter_gaussian(tmag, old_params), /ylog, xrange = [-23, -19], /xsty
;oplot, tmag2, schechter_mag(tmag2, blan), color = !red
;oplot, tmag2, tphi*phi_rescale, color = !blue

end
