pro range, x, trunc=trunc, sig=sig
print, min(x), max(x)
if keyword_set(trunc) then print, medscat(x,0.01), medscat(x, 0.99)
if keyword_set(sig) then print, medscat(x,0.167), medscat(x, 0.833)
end
