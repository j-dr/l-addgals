function get_bcc_truth_structure

obs_fullone=create_struct('id',0LL,'index', 0L, 'ecatid', 0L, 'coeffs', fltarr(5), $
                          'tmag',fltarr(5),'omag',fltarr(5), $
                          'flux', fltarr(5), 'ivar', fltarr(5), $
                          'omagerr',fltarr(5),'amag', fltarr(5), $
                          'ra',0.0,'dec',0.0,'z',0.0, $
                          'haloid',0L,'rhalo',0.0, 'm200', 0., 'ngals', 0L, $
                          'r200',0.0,'central',0, 'tra', 0., 'tdec', 0., $
                          'epsilon', fltarr(2), 'gamma1', 0., 'gamma2', 0., $
                          'kappa', 0., 'mu', 0., 'lmag', fltarr(5), 'mag_u', 0.,$
                          'te', fltarr(2), 'tsize', 0., 'size', 0., 'size_subaru', 0., $
                          'ArborZ', 0., $
                          'ArborZ_err', 0., 'ANNZ', 0., 'ANNZ_err', 0.,$
                          'photoz_Gaussian', 0., $
                          'px', 0., 'py', 0., 'pz', 0., $
                          'vx', 0., 'vy', 0., 'vz', 0., $
			  'pstar', 0., 'pqso', 0., $
			  'n_sersic', 0., 'r50', 0.)

return, obs_fullone

end
