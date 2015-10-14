function exponential_power_distribution_cdf, x, alpha, beta

cdf = igamma(1/beta,(abs(x)/alpha)^beta)/(2*gamma(1./beta))
ii = where(x lt 0, count)
if (count gt 0) then cdf[ii] *= -1
cdf += 0.5

return, cdf

end


