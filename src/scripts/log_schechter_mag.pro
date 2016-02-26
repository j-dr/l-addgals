FUNCTION log_schechter_mag, mag, p
;Mstar = p[0]
;phistar = p[1]
;alpha = p[2]
return, alog10(0.4*alog(10.0)*p[0]*10^(-0.4*(mag-p[2])*(p[1]+1))*exp(-1*10^(-0.4*(mag-p[2]))))
end
