function get_bcc_obs_structure

 obs_visibleone=create_struct('id',0LL, 'index', 0L, 'mag_u', 0.,$
                               'mag_g',0.,'mag_r',0.,'mag_i',0.,'mag_z',0.,$
                               'mag_y',0., 'magerr_g', 0., 'magerr_r', 0., $
                               'magerr_i', 0., 'magerr_z', 0., 'magerr_y', 0., $
                               'flux_g',0.,'flux_r',0.,'flux_i',0.,'flux_z',0.,'flux_y',0.,$
                               'ivar_g', 0., 'ivar_r', 0., 'ivar_i', 0., $
                               'ivar_z',0., 'ivar_y', 0., $
                               'ra',0.0,'dec',0.0,$
			       'epsilon1', 0., 'epsilon2', 0., 'size', 0., $
			       'pstar', 0.0, 'pqso', 0.0, $
                               'ArborZ', 0., $
                               'ArborZ_err', 0., 'ANNZ', 0., 'ANNZ_err', 0.,$
                               'photoz_Gaussian', 0.)

return, obs_visibleone
end
