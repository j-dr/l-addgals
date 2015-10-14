.comp mask_bcc_pixel.pro

resolve_all, skip_routines=['trnlog', 'dellog', 'red']

save, filename='mask_bcc_pixel.sav', /routines

exit

