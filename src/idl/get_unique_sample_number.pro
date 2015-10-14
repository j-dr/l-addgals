function get_unique_sample_number, ind, num

nind = N_ELEMENTS(ind)

irand = randomu(seed, nind)
srt = sort(irand)
jj = ind[srt[0:num-1]]

return, jj

end

