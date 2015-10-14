.comp make_des_220_sq_deg_catalog.pro

resolve_all, skip_routines=['trnlog', 'dellog', 'red']

save, filename='make_des_220_sq_deg_catalog.sav', /routines

exit

