function lasdamas_z, num, verbose=verbose

zred = 1./(0.075*1.0265097^float(num))-1

if KEYWORD_SET(verbose) then begin
   print, "zred = 0.075*1.0265097^<num>"
endif

return, zred

end
