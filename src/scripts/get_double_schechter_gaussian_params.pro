function get_double_schechter_gaussian_params, err=err

mstar = -19.88
mstar_err = 0.03
phi1 = 0.0156
phi1_err = 0.0005
alpha1 = -0.166
alpha1_err = 0.041
phi2 = 0.00671
phi2_err = 0.00029
alpha2 = -1.523
alpha2_err = 0.01
phi3 = 3.08e-5
phi3_err = 3.24e-5
Mhi = -21.72
Mhi_err = 0.52
sigma_hi = 0.484
sigma_hi_err = 0.192

params = [phi1, alpha1, phi2, alpha2, mstar, phi3, mhi, sigma_hi]
err = [phi1_err, alpha1_err, phi2_err, alpha2_err, mstar_err, $
	phi3_err, mhi_err, sigma_hi_err]

return, params

end
