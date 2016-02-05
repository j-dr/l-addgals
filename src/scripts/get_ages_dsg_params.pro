;;;short routine to return best-fitting parameters for the AGES LF

function get_ages_dsg_params, zred, remove_passive_evolution=remove_passive_evolution, Q = Q

;if not KEYWORD_SET(Q) then Q = -1.49
;if not KEYWORD_SET(Q) then Q = 3.16
if not KEYWORD_SET(Q) then Q = 3.02428

;;;our AGES parameters in redshift bins
z_ages = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65]
n_ages = N_ELEMENTS(z_ages)
mstar_ages = [-20.27, -20.44, -20.81, -20.81, -20.99, -21.29, -21.38]
mstar_ages_err = [0.0, 0.05, 0.04, 0.03, 0.04, 0.08, 0.06]
phistar_ages = [1.60, 1.49, 1.52, 1.24, 1.44, 1.08, 1.05]*1e-2;*0.97366
phistar_ages_err_percent = [0.0, 0.011, 0.05, 0.07, 0.07, 0.1, 0.14]
phistar_ages_err = phistar_ages*phistar_ages_err_percent
alpha = -1.05
;mstar0 = -20.27
;phistar0 = 1.65e-2 ;;not actually used anywhere
mstar0 = -20.310 ;;value from Rachel

;;;do we just interpolate between tabulated ages parameters?
if (zred le z_ages(n_ages-1)) then begin
  ka = min(abs(zred - z_ages), loc)
  if (loc eq n_ages - 1) then loc--
  slope_mstar = (mstar_ages[loc+1] - mstar_ages[loc]) / $
		(z_ages[loc+1] - z_ages[loc])
  mstar = mstar_ages[loc] + slope_mstar*(zred-z_ages[loc]) 
  slope_phistar = alog10(phistar_ages[loc+1]/phistar_ages[loc]) / $
		  (z_ages[loc+1] = z_ages[loc])
  phistar = 10.^(alog10(phistar_ages[loc]) + slope_phistar*(zred-z_ages[loc]))
endif else begin ;;here we extrapolate our model beyond the largest redshift
  phistar = 0.0168*10.^(0.4*(-0.972)*zred)
  mstar = -20.44 - Q*(1./(1+zred)-1/1.1)
endelse

;;;redo the parameters using our smooth model
phistar = 10.^(-1.79574 + (-0.266409*zred))
mstar = -23.1907 + Q*(1./(1+zred))


if (KEYWORD_SET(remove_passive_evolution)) then begin
	mstar -= Q*(1./(1+zred)-1./1.1)
;	mstar0 -= Q*(zred-0.1) ;;don't worry about this -- it comes out in the shift
	if (zred ge 0.0) then mstar = -20.44
;	phistar = 0.0168*10.^(0.4*(-0.972)*zred)
endif

params = get_double_schechter_gaussian_params()
phistar_rat = phistar/params[0]
mr_shift = mstar - mstar0
params[0] *= phistar_rat
params[2] *= phistar_rat
params[5] *= phistar_rat
params[4] += mr_shift
params[6] += mr_shift

return, params

end


