pro calc_elipticities, mag, e1, e2, emagmin=emagmin, emagmax=emagmax

;These parameters define Joerg's fit to the shape distribution
if not KEYWORD_SET(emagmin) then emagmin = 21.38
if not KEYWORD_SET(emagmax) then emagmax = 26.81
if not KEYWORD_SET(erefmag) then erefmag = 24.0
if not KEYWORD_SET(epsmax) then epsmax = 1.0
if not KEYWORD_SET(p_0) then p_0 = 1.2565
if not KEYWORD_SET(p_1) then p_1 =  0.0937
if not KEYWORD_SET(p_2) then p_2 =  -0.0049
if not KEYWORD_SET(p_3) then p_3 = -0.0029
if not KEYWORD_SET(sigma_0) then sigma_0 = 0.2679
if not KEYWORD_SET(sigma_1) then sigma_1 = 0.0415
if not KEYWORD_SET(sigma_2) then sigma_2 = 0.0010
if not KEYWORD_SET(sigma_3) then sigma_3 = -0.0010

epsmax2 = epsmax*epsmax
sqrt2 = sqrt(2.)
ng = N_ELEMENTS(mag)
mymag = mag

ii = where(mymag lt emagmin, count) 
if (count gt 0) then mymag[ii] = emagmin
ii = where(mymag gt emagmax) 
if (count gt 0) then mymag[ii] = emagmax

mymag -= erefmag

beta = p_0 + mymag*(p_1 + mymag*(p_2 + mymag*p_3))
sigma_p = sigma_0 + mymag*(sigma_1 + mymag*(sigma_2 + mymag*sigma_3))
alpha = sigma_p*beta^(1./beta)


nrand = 1000LL
xval = (0.5+findgen(nrand))/(0.5*nrand)*epsmax-epsmax
ran = randomu(seed,ng)
eps = fltarr(ng)
for i = 0L, ng - 1 do begin
  cdf = exponential_power_distribution_cdf(xval, (1.14/0.85)*sqrt(2)*alpha[i], beta[i])
  ka = min(abs(cdf-ran[i]), loc)
  eps[i] = xval[loc]
endfor


asdf

phi = randomu(seed,ng)*2*!PI - !PI
e1 = e*cos(2*phi)
e2 = e*sin(2*phi)

end
