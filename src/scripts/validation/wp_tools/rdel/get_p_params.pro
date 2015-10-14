function get_p_params, paramfile

; read everythign form the paramfile
readcol, paramfile, label, value, format = 'a,f'

ind1 = where(label eq 'p0')
ind2 = where(label eq 'p1')
ind3 = where(label eq 'p2')
ind4 = where(label eq 'p3')
ind5 = where(label eq 'pz1')
ind6 = where(label eq 'pz2')
ind7 = where(label eq 'pz3')

params = value([ind1, ind2, ind3, ind4, ind5, ind6, ind7])

return, params

end
