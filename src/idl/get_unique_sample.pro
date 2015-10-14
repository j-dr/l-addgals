function get_unique_sample, ind, frac

nind = N_ELEMENTS(ind)

irand = randomu(seed, nind)
jj = where(irand lt frac)

return, jj

end
