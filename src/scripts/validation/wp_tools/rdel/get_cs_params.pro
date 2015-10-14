function get_cs_params, paramfile

; read everythign form the paramfile
readcol, paramfile, label, value, format = 'a,f'

ind1 = where(label eq 'cs0')
ind2 = where(label eq 'cs1')
ind3 = where(label eq 'cs2')
ind4 = where(label eq 'cs3')
ind5 = where(label eq 'csz1')
ind6 = where(label eq 'csz2')

params = value([ind1, ind2, ind3, ind4, ind5, ind6])

return, params

end
