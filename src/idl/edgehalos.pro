pro edgehalos, h
maxz=0.34
;maxz = max(g.z)
print, maxz
edge = where(h.halopx lt 0+h.r200 or h.halopy lt 0+h.r200 or h.halopz lt 0+h.r200 or h.z gt maxz-0.00333, comp=comp, count)
print, count, n_elements(h)
if (count gt 0) then h[edge].edge = 1
end
