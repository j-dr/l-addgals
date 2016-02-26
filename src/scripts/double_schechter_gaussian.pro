function double_schechter_gaussian, mr, params

phi1 = params[0]
alpha1 = params[1]
phi2 = params[2]
alpha2 = params[3]
mstar = params[4]
phi3 = params[5]
Mhi = params[6]
sigma_hi = params[7]

phi = 0.4*alog(10.)*exp(-10.^(-0.4*(mr-mstar))) * ($
	phi1*10.^(-0.4*(mr-mstar)*(alpha1+1)) + $
	phi2*10.^(-0.4*(mr-mstar)*(alpha2+1))) + $
      	phi3/sqrt(2*!pi*sigma_hi^2.)*exp(-(mr-Mhi)^2./(2*sigma_hi^2))

return, phi

end
