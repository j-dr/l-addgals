function get_cm_params, paramfile

; read everythign form the paramfile
readcol, paramfile, label, value, format = 'a,f'

ind1 = where(label eq 'cm0')
ind2 = where(label eq 'cm1')
ind3 = where(label eq 'cm2')
ind4 = where(label eq 'cm3')
ind5 = where(label eq 'cm4')
ind6 = where(label eq 'cmz1')
ind7 = where(label eq 'cmz2')
ind8 = where(label eq 'cmz3')

params = value([ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8])

return, params

end
