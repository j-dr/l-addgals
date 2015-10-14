.comp finalize_bcc_catalog.pro

resolve_all, skip_routines=['trnlog', 'dellog', 'red']

save, filename='finalize_bcc_catalog.sav', /routines

exit
