;program: split_des_file
;purpose: cuts one of the observed_Xyear catalog files into two files
;         with equal area
;
;parameters:
;         base: The name of the file to split, i.e. 'hv_POv1.06b_oct'
;          num: The number of the file, i.e. 1 for file 'hv_POv1.06b_oct.01_observed_5year'
; new_base_num: The first of the two numbers for the split file -- the
;               two files will be numbered new_base_num and
;               new_base_num+1.  Make sure you don't accidentally
;               write over the old file.  

pro split_des_file, base, num, new_base_num

if (num eq new_base_num) then begin
  print, " You don't really want to write over the file, do you?"
  print, "   Specify nem != new_base_num"
  stop
endif

num_sfx = '0'+strcompress(string(num), /remove_all)
num_sfx = strmid(num_sfx, strlen(num_sfx)-2,2)
base += '.'
dim_sfx = ['', '_dim']
gal_sfx = ['_observed_1year','_observed_5year']

for gi = 0,1 do begin
for dim = 0, 1 do begin
  g = mrdfits(base+num_sfx+gal_sfx(gi)+dim_sfx(dim)+'.fit', 1)

  ;;;determine where we split
  min_ra = min(g.ra)
  max_ra = max(g.ra)
  mid_ra = 0.5*(min_ra+max_ra)

  ;;;split the galaxies and halos
  g_ind1 = where(g.ra lt mid_ra)
  g_ind2 = where(g.ra ge mid_ra)

  ;;;write out the information
  outnum1 = '0'+strcompress(string(new_base_num), /remove_all)
  outnum1 = strmid(outnum1, strlen(outnum1)-2,2)
  outnum2 = '0'+strcompress(string(new_base_num+1), /remove_all)
  outnum2 = strmid(outnum2, strlen(outnum2)-2,2)

  mwrfits, g(g_ind1), base+outnum1+gal_sfx(gi)+dim_sfx(dim)+'.fit', /create
  mwrfits, g(g_ind2), base+outnum2+gal_sfx(gi)+dim_sfx(dim)+'.fit', /create

  g = 0
endfor
endfor


end
 
