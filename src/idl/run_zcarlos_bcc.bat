.comp run_zcarlos_bcc.pro

resolve_all, skip_routines=['trnlog', 'dellog', 'red']

save, filename='run_zcarlos_bcc.sav', /routines

exit

