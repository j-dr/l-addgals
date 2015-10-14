.comp run_zcarlos_bossanova.pro

resolve_all, skip_routines=['trnlog', 'dellog', 'red']

save, filename='run_zcarlos_bossanova.sav', /routines

exit

