pro get_dr8_full_sky_training_set, path, outpath=outpath, $
	truth_base=truth_base, sdss_base = sdss_base, $
	outfile1=outfile1, outfile2=outfile2

if not KEYWORD_SET(truth_base) then truth_base = 'PO_Aardvark_truth'
if not KEYWORD_SET(sdss_base) then sdss_base = 'PO_Aardvark_sdss_mag'
if not KEYWORD_SET(outpath) then outpath = path
if not KEYWORD_SET(outfile1) then outfile1 = outpath+'/PO_Aardvark_DR8_training_set.fit'
if not KEYWORD_SET(outfile2) then outfile2 = outpath+'/PO_Aardvark_DR8_training_set_sdss_mag.fit'
frac_zcosmos = 0.18
frac_cnoc2 = 0.18
frac_deep2 = 0.43
frac_vvds = 0.053
frac_primus_cosmos = 0.48
frac_primus_deep2_1 = 0.48
frac_primus_deep2_2 = 0.48
frac_primus_xmm = 0.48
frac_2slaq = 0.5

;;;the zcosmos sample
print, 'Getting zCOSMOS sample...'
file = path+'/truth/'+truth_base+'.124.fit'
file2 = path+'/DR8/'+sdss_base+'.124.fit'
g = mrdfits(file, 1)
s = mrdfits(file2, 1)
sp = cut_to_dr8_training_limits(s, ii=ii)
gp = g[ii]
ii_zcosmos =        where(gp.ra gt 140 and gp.ra lt 141  and gp.dec gt 41   and gp.dec lt 43.3, n_zcosmos)
jj = get_unique_sample(ii_zcosmos, frac_zcosmos)
ii_zcosmos = ii_zcosmos[jj]
g_train = gp[ii_zcosmos]
s_train = sp[ii_zcosmos]

;;;the CNOC2 sample
print, 'Getting CNOC2 sample...'
file = path+'/truth/'+truth_base+'.181.fit'
file2 = path+'/DR8/'+sdss_base+'.181.fit'
g = mrdfits(file, 1)
s = mrdfits(file2, 1)
sp = cut_to_dr8_training_limits(s, ii=ii)
gp = g[ii]
ii_cnoc2 =          where(gp.ra gt 61 and gp.ra lt 62.145    and gp.dec gt 29   and gp.dec lt 29.4, n_conc2)
jj = get_unique_sample(ii_cnoc2, frac_cnoc2)
ii_cnoc2 = ii_cnoc2[jj]
g_train = [g_train, gp[ii_cnoc2]]
s_train = [s_train, sp[ii_cnoc2]]

;;;the DEEP2 sample
print, 'Getting DEEP2 sample...'
file = path+'/truth/'+truth_base+'.187.fit'
file2 = path+'/DR8/'+sdss_base+'.187.fit'
g = mrdfits(file, 1)
s = mrdfits(file2, 1)
sp = cut_to_dr8_training_limits(s, ii=ii)
gp = g[ii]
ii_deep2 =          where(gp.ra gt 129 and gp.ra lt 130.14    and gp.dec gt 29   and gp.dec lt 29.4, n_deep2)
jj = get_unique_sample(ii_deep2, frac_deep2)
ii_deep2 = ii_deep2[jj]
g_train = [g_train, gp[ii_deep2]]
s_train = [s_train, sp[ii_deep2]]

;;;the VVDS sample
print, 'Getting VVDS sample...'
file = path+'/truth/'+truth_base+'.209.fit'
file2 = path+'/DR8/'+sdss_base+'.209.fit'
g = mrdfits(file, 1)
s = mrdfits(file2, 1)
sp = cut_to_dr8_training_limits(s, ii=ii)
gp = g[ii]
ii_vvds =           where(gp.ra gt 10 and gp.ra lt 12.2    and gp.dec gt 23   and gp.dec lt 25, n_vvds)
jj = get_unique_sample(ii_vvds, frac_vvds)
ii_vvds = ii_vvds[jj]
g_train = [g_train, gp[ii_vvds]]
s_train = [s_train, sp[ii_vvds]]

;;;the main PRIMUS samples
print, 'Getting PRIMUS samples...'
file = path+'/truth/'+truth_base+'.245.fit'
file2 = path+'/DR8/'+sdss_base+'.245.fit'
g = mrdfits(file, 1)
s = mrdfits(file2, 1)
sp = cut_to_dr8_training_limits(s, ii=ii)
gp = g[ii]
ii_primus_cosmos =  where(gp.ra gt 62 and gp.ra lt 63 and gp.dec gt 18   and gp.dec lt 19.16, n_primus_c)
jj = get_unique_sample(ii_primus_cosmos, frac_primus_cosmos)
ii_primus_cosmos = ii_primus_cosmos[jj]
g_train = [g_train, gp[ii_primus_cosmos]]
s_train = [s_train, sp[ii_primus_cosmos]]

;;;PRIMUS field 1
file = path+'/truth/'+truth_base+'.275.fit'
file2 = path+'/DR8/'+sdss_base+'.275.fit'
g = mrdfits(file, 1)
s = mrdfits(file2, 1)
sp = cut_to_dr8_training_limits(s, ii=ii)
gp = g[ii]
ii_primus_deep2_1 = where(gp.ra gt 33 and gp.ra lt 34 and gp.dec gt 14   and gp.dec lt 14.72, n_primus_d1)
jj = get_unique_sample(ii_primus_deep2_1, frac_primus_deep2_1)
ii_primus_deep2_1 = ii_primus_deep2_1[jj]
g_train = [g_train, gp[ii_primus_deep2_1]]
s_train = [s_train, sp[ii_primus_deep2_1]]

;;;PRIMUS field 2
file = path+'/truth/'+truth_base+'.317.fit'
file2 = path+'/DR8/'+sdss_base+'.317.fit'
g = mrdfits(file, 1)
s = mrdfits(file2, 1)
sp = cut_to_dr8_training_limits(s, ii=ii)
gp = g[ii]
ii_primus_deep2_2 = where(gp.ra gt 151 and gp.ra lt 152    and gp.dec gt 9 and gp.dec lt 9.71, n_primus_d2)
jj = get_unique_sample(ii_primus_deep2_2, frac_primus_deep2_2)
ii_primus_deep2_2 = ii_primus_deep2_2[jj]
g_train = [g_train, gp[ii_primus_deep2_2]]
s_train = [s_train, sp[ii_primus_deep2_2]]

;;;PRIMUS XMM
file = path+'/truth/'+truth_base+'.350.fit'
file2 = path+'/DR8/'+sdss_base+'.350.fit'
;file = path+'/truth/'+truth_base+'.219.fit'
;file2 = path+'/DR8/'+sdss_base+'.219.fit'
g = mrdfits(file, 1)
s = mrdfits(file2, 1)
sp = cut_to_dr8_training_limits(s, ii=ii)
gp = g[ii]
ii_primus_xmm =     where(gp.ra gt 157 and gp.ra lt 158.71    and gp.dec gt 4   and gp.dec lt 6, n_primus_xmm)
;ii_primus_xmm =     where(gp.ra gt 123 and gp.ra lt 125    and gp.dec gt 24   and gp.dec lt 25.59, n_primus_xmm)
jj = get_unique_sample(ii_primus_xmm, frac_primus_xmm)
ii_primus_xmm = ii_primus_xmm[jj]
g_train = [g_train, gp[ii_primus_xmm]]
s_train = [s_train, sp[ii_primus_xmm]]

;;;SDSS DR5
print, 'Getting DR5 sample...'
target_dr5 = 100000L
n_dr5 = 0L
for i = 0, 400 do begin
  file = path+'/truth/'+truth_base+'.'+strcompress(string(i),/remove_all)+'.fit'
  file2 = path+'/DR8/'+sdss_base+'.'+strcompress(string(i),/remove_all)+'.fit'
  if (file_test(file) eq 0) then continue
  g = mrdfits(file,1)
  s = mrdfits(file2, 1)
  sp = cut_to_dr8_training_limits(s, ii=ii)
  gp = g[ii]
  ii_dr5 = where(sp.omag[2] lt 17.8, tn_dr5)
  if (n_dr5 + tn_dr5 gt target_dr5) then begin
    n_to_keep = target_dr5 - n_dr5
    ii_dr5 = ii_dr5[0:n_to_keep-1]
  endif
  n_dr5 += tn_dr5
  g_train = [g_train,gp[ii_dr5]]
  s_train = [s_train,sp[ii_dr5]]
  if (n_dr5 gt target_dr5) then break
endfor
istart = i+1

print, 'Getting LRG sample...'
target_lrg = 20000L
n_lrg = 0L
for i = istart, 400 do begin
  file = path+'/truth/'+truth_base+'.'+strcompress(string(i),/remove_all)+'.fit'
  file2 = path+'/DR8/'+sdss_base+'.'+strcompress(string(i),/remove_all)+'.fit'
  if (file_test(file) eq 0) then continue
  g = mrdfits(file,1)
  s = mrdfits(file2, 1)
  sp = cut_to_dr8_training_limits(s, ii=ii)
  gp = g[ii]
  gmr = sp.omag[1] - sp.omag[2]
  rmi = sp.omag[2] - sp.omag[3]
  ii_lrg = where(sp.omag[2] lt 19.5 and gmr lt 1.3+0.25*(rmi), tn_lrg)
  if (n_lrg + tn_lrg gt target_lrg) then begin
    n_to_keep = target_lrg - n_lrg
    ii_lrg = ii_lrg[0:n_to_keep-1]
  endif
  g_train = [g_train, gp[ii_lrg]]
  s_train = [s_train, sp[ii_lrg]]
  n_lrg += tn_lrg
  if (n_lrg ge target_lrg) then break
endfor
istart = i+1

print, 'Getting 2SLAQ sample...'
target_2slaq = 8633
n_2slaq = 0L
for i = istart, 400 do begin
  file = path+'/truth/'+truth_base+'.'+strcompress(string(i),/remove_all)+'.fit'
  file2 = path+'/DR8/'+sdss_base+'.'+strcompress(string(i),/remove_all)+'.fit'
  if (file_test(file) eq 0) then continue
  g = mrdfits(file,1)
  s = mrdfits(file2, 1)
  sp = cut_to_dr8_training_limits(s, ii=ii)
  gp = g[ii]
  gmr = sp.omag[1] - sp.omag[2]
  rmi = sp.omag[2] - sp.omag[3]
  ii_2slaq = where(rmi-gmr/8 ge 0.55 and $
		 0.7*gmr + 1.2*(rmi-0.18) ge 1.6 and $
		 16.5 le sp.omag[3] and sp.omag[3] lt 19.6 and $
		 0.5 lt gmr and gmr lt 3 and $
		 rmi lt 2, tn_2slaq)
  jj = get_unique_sample(ii_2slaq, 0.5)
  ii_2slaq = ii_2slaq[jj]
  tn_2slaq = N_ELEMENTS(ii_2slaq)
  if (n_2slaq + tn_2slaq gt target_2slaq) then begin
    n_to_keep = target_2slaq - n_2slaq
    ii_2slaq = ii_2slaq[0:n_to_keep-1]
  endif
  g_train = [g_train, gp[ii_2slaq]]
  s_train = [s_train, sp[ii_2slaq]]
  n_2slaq += tn_2slaq
  if (n_2slaq ge target_2slaq) then break
endfor

add_tags, s_train, ['ra', 'dec', 'z'], ['0.', '0.', '0.'], s_out
s_out.z = g_train.z
s_out.ra = g_train.ra
s_out.dec = g_train.dec
mwrfits, g_train, outfile1, /create
mwrfits, s_out, outfile2, /create

end
