.comp run_zcarlos_dr8.pro

resolve_all, skip_routines=['trnlog', 'dellog', 'red']

save, filename='run_zcarlos_dr8.sav', /routines

exit

