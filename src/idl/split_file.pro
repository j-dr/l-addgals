pro split_file, base, num, new_base_num

num_sfx = '0'+strcompress(string(num), /remove_all)
num_sfx = strmid(num_sfx, strlen(num_sfx)-2,2)
base += '.'
dim_sfx = ['', '_dim']
gal_sfx = '_galaxies'
halo_sfx = '_halos'
survey_sfx = ['_des', '_vista', '_deep', '_johnson', '_ubv_deep']

for dim = 1, 1 do begin
  g = mrdfits(base+num_sfx+gal_sfx+dim_sfx(dim)+'.fit', 1)
  h = mrdfits(base+num_sfx+halo_sfx+dim_sfx(dim)+'.fit', 1)

  ;;;determine where we split
  min_ra = min(g.ra)
  max_ra = max(g.ra)
  mid_ra = 0.5*(min_ra+max_ra)

  ;;;rename the nonedge tag to edge
  g = rename_tags(g, 'nonedge', 'edge')
  des = mrdfits(base+num_sfx+survey_sfx(0)+dim_sfx(dim)+'.fit', 1)
  if (dim eq 1) then begin
    ii = where(des.omag(2) le 24)
    g = g[ii]
    des = des[ii]
  endif

  ;;;split the galaxies and halos
  g_ind1 = where(g.ra lt mid_ra)
  g_ind2 = where(g.ra ge mid_ra)
  h_ind1 = where(h.ra lt mid_ra)
  h_ind2 = where(h.ra ge mid_ra)

  ;;;find the new edge galaxies
  ;;;start by sorting the halo ids for easy lookup
  HaloArray = lonarr(max(h.haloid)+1)
  for hi = 0L, N_ELEMENTS(h) - 1 do begin
    HaloArray(h(hi).haloid) = hi
  endfor

  ;;;update the edge tags
  for i = 0L, N_ELEMENTS(g)-1 do begin
    gal_file = 0
    if (g(i).ra ge mid_ra) then gal_file = 1
    halo_file = 0
    if (h(HaloArray(g(i).haloid)).ra ge mid_ra) then halo_file = 1
    if (gal_file ne halo_file) then begin
      g(i).edge = halo_file+new_base_num
      h(HaloArray(g(i).haloid)).edge = 1
    endif else begin
      g(i).edge = 0
      h(HaloArray(g(i).haloid)).edge = 0
    endelse
  endfor

  ;;;write out the information
  outnum1 = '0'+strcompress(string(new_base_num), /remove_all)
  outnum1 = strmid(outnum1, strlen(outnum1)-2,2)
  outnum2 = '0'+strcompress(string(new_base_num+1), /remove_all)
  outnum2 = strmid(outnum2, strlen(outnum2)-2,2)

  mwrfits, g(g_ind1), base+outnum1+gal_sfx+dim_sfx(dim)+'.fit', /create
  mwrfits, g(g_ind2), base+outnum2+gal_sfx+dim_sfx(dim)+'.fit', /create
  mwrfits, h(h_ind1), base+outnum1+halo_sfx+dim_sfx(dim)+'.fit', /create
  mwrfits, h(h_ind2), base+outnum2+halo_sfx+dim_sfx(dim)+'.fit', /create
  mwrfits, des(g_ind1), base+outnum1+survey_sfx(0)+dim_sfx(dim)+'.fit', /create
  mwrfits, des(g_ind2), base+outnum2+survey_sfx(0)+dim_sfx(dim)+'.fit', /create

  ;;;free up space because we'll be readin in the other files
  g = 0
  h = 0
  des = 0
  ;;;start at 1 b/c we already did des for the magnitdue cut
  for survey = 1, 3 do begin
    g_survey = mrdfits(base+num_sfx+survey_sfx(survey)+dim_sfx(dim)+'.fit', 1)
    ;;;if dim we enforce a magnitude cut
    if (dim eq 1) then begin
      g_survey = g_survey[ii]
    endif
    mwrfits, g_survey(g_ind1), base+outnum1+survey_sfx(survey)+dim_sfx(dim)+'.fit', /create
    mwrfits, g_survey(g_ind2), base+outnum2+survey_sfx(survey)+dim_sfx(dim)+'.fit', /create
    g_survey = 0
  endfor

endfor

end
 
