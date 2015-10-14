pro cal_sizes, mag, size, xi_0=xi_0, xi_1=xi_1, xi_2=xi_2, xi_3=xi_3, xi_4=xi_4,$
	alpha_0=alpha_0,alpha_1=alpha_1,alpha_2=alpha_2,alpha_3=alpha_3,alpha_4=alpha_4, $
	kappa_0=kappa_0,kappa_1=kappa_1,kappa_2=kappa_2,kappa_3=kappa_3,kappa_4=kappa_4

if not KEYWORD_SET(xi_0) then xi_0 = 
if not KEYWORD_SET(xi_1) then xi_1 = 
if not KEYWORD_SET(xi_2) then xi_2 = 
if not KEYWORD_SET(xi_3) then xi_3 = 
if not KEYWORD_SET(xi_4) then xi_4 = 
if not KEYWORD_SET(alpha_0) then alpha_0 = 
if not KEYWORD_SET(alpha_1) then alpha_1 = 
if not KEYWORD_SET(alpha_2) then alpha_2 = 
if not KEYWORD_SET(alpha_3) then alpha_3 = 
if not KEYWORD_SET(alpha_4) then alpha_4 = 
if not KEYWORD_SET(kappa_0) then kappa_0 = 
if not KEYWORD_SET(kappa_1) then kappa_1 = 
if not KEYWORD_SET(kappa_2) then kappa_2 = 
if not KEYWORD_SET(kappa_3) then kappa_3 = 
if not KEYWORD_SET(kappa_4) then kappa_4 = 

xi = xi_0 + mag*(xi_1 + mag*(xi_2 + mag*(xi_3 + mag*xi_4)))
alpha = alpha_0 + mag*(alpha_1 + mag*(alpha_2 + mag*(alpha_3 + mag*alpha_4)))
kappa = kappa_0 + mag*(kappa_1 + mag*(kappa_2 + mag*(kappa_3 + mag*kappa_4)))


