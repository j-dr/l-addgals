function get_fs_params, paramfile

; read everythign form the paramfile
readcol, paramfile, label, value, format = 'a,f'

ind1 = where(label eq 'fs0')
ind2 = where(label eq 'fs1')
ind3 = where(label eq 'fs2')
ind4 = where(label eq 'fs3')
ind5 = where(label eq 'fs4')
ind6 = where(label eq 'fsz1')
ind7 = where(label eq 'fsz2')
ind8 = where(label eq 'fsz3')

params = value([ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8])

return, params

end
