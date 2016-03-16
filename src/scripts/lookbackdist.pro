;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION: LookbackDist
; PURPOSE:  Returns lookback dist in Mpc to a specific a given 
;           OmegaM and OmegaL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function LookbackDist, z, OmegaM, OmegaL, h = h

a = 1./(1+z)
c = 2.99792458e5 ;speek of light in km/s
da = 1e-6
r = 0.
Thisa = 1.
i = 0.
if (ARG_PRESENT(h) eq 0) then h = 1.0

while (Thisa gt a) do begin
  Thisa = 1.-da*i
  ThisH = 100*h*sqrt(OmegaM/Thisa^3 + OmegaL)
  r += 1./(ThisH*Thisa*Thisa)
  i++
endwhile

r *= da*c

return, r

end
