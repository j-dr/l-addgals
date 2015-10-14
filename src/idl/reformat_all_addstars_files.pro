pro reformat_all_addstars_files

;;;out input/output directories
file_list = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Stars/unique_files.txt'
outpath = '/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Stars/temp_files/'

;;;read a list of the individual catalog files
maxfiles = 10000
nfiles = 0L
infile = strarr(maxfiles)
entry = 'asdf'
openr,1,file_list
while (eof(1) eq 0) do begin
  readf,1,entry
  infile[nfiles] = entry
  nfiles++
  if (nfiles eq maxfiles) then begin
    print, 'Error!  need to increase maxfiles!'
    break
  endif
endwhile
close,1

;;;loop through all the individual addstars files
print, 'Now reformatting '+strcompress(string(nfiles))+' files...'
for i = 0L, nfiles - 1 do begin
  print, 'processing file '+strcompress(string(i))+': '+infile[i]
  reformat_addstars_file, infile[i], outpath
endfor

end
