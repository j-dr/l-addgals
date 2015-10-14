function get_fm_params, paramfile

; read everythign form the paramfile
readcol, paramfile, label, value, format = 'a,f'

ind1 = where(label eq 'fm0')
ind2 = where(label eq 'fm1')
ind3 = where(label eq 'fm2')
ind4 = where(label eq 'fm3')
ind5 = where(label eq 'fmz1')
ind6 = where(label eq 'fmz2')

params = value([ind1, ind2, ind3, ind4, ind5, ind6])

return, params

end
