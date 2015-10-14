function get_dr8_training_set, g

ii1 = where(g.omag[2] lt 21.8, count1)

for i = 22.1, 32.1, 0.1 do begin
  ii2 = where(g[ii1].omag[0] lt i and g[ii1].omag[1] lt 29 and $
	      g[ii1].omag[3] lt 29 and g[ii1].omag[4] lt 29, count2)
  if (float(count2)/float(count1) gt 0.8) then break
endfor

ii = ii1[ii2]
gp = g[ii]

ii_zcosmos =        where(gp.ra gt 12 and gp.ra lt 13.5  and gp.dec gt 37   and gp.dec lt 38.44, n_zcosmos)
ii_cnoc2 =          where(gp.ra gt 28 and gp.ra lt 29    and gp.dec gt 38   and gp.dec lt 38.51, n_conc2)
ii_deep2 =          where(gp.ra gt 14 and gp.ra lt 15    and gp.dec gt 45   and gp.dec lt 45.57, n_deep2)
ii_vvds =           where(gp.ra gt 25 and gp.ra lt 27    and gp.dec gt 38   and gp.dec lt 40.59, n_vvds)
ii_primus_cosmos =  where(gp.ra gt 21 and gp.ra lt 22.45 and gp.dec gt 40   and gp.dec lt 41, n_primus_c)
ii_primus_deep2_1 = where(gp.ra gt 23 and gp.ra lt 24.09 and gp.dec gt 44   and gp.dec lt 44.9, n_primus_d1)
ii_primus_deep2_2 = where(gp.ra gt 17 and gp.ra lt 18    and gp.dec gt 38.5 and gp.dec lt 39.4, n_primus_d2)
ii_primus_xmm =     where(gp.ra gt 20 and gp.ra lt 22    and gp.dec gt 46   and gp.dec lt 48.5, n_primus_xmm)
gmr = gp.omag[1] - gp.omag[2]
rmi = gp.omag[2] - gp.omag[3]
ii_dr5 = where(gp.omag[2] lt 17.8, n_dr5)
ii_lrg = where(gp.omag[2] lt 19.5 and gmr lt 1.3+0.25*(rmi), n_lrg)
ii_2slaq = where(rmi-gmr/8 ge 0.55 and $
		 0.7*gmr + 1.2*(rmi-0.18) ge 1.6 and $
		 16.5 le gp.omag[3] and gp.omag[3] lt 19.6 and $
		 0.5 lt gmr and gmr lt 3 and $
		 rmi lt 2 and $
		 gp.ra gt 10 and gp.ra lt 30 and $
		 gp.dec gt 35 and gp.dec lt 47, n_2slaq)

frac_zcosmos = 0.18
frac_cnoc2 = 0.18
frac_deep2 = 0.43
frac_vvds = 0.053
frac_primus_cosmos = 0.48
frac_primus_deep2_1 = 0.48
frac_primus_deep2_2 = 0.48
frac_primus_xmm = 0.48
frac_dr5 = 100000./n_dr5
frac_lrg = 20000./n_lrg
frac_2slaq = 0.5

jj = get_unique_sample(ii_zcosmos, frac_zcosmos)
ii_zcosmos = ii_zcosmos[jj]
jj = get_unique_sample(ii_cnoc2, frac_cnoc2)
ii_cnoc2 = ii_cnoc2[jj]
jj = get_unique_sample(ii_deep2, frac_deep2)
ii_deep2 = ii_deep2[jj]
jj = get_unique_sample(ii_vvds, frac_vvds)
ii_vvds = ii_vvds[jj]
jj = get_unique_sample(ii_primus_cosmos, frac_primus_cosmos)
ii_primus_cosmos = ii_primus_cosmos[jj]
jj = get_unique_sample(ii_primus_deep2_1, frac_primus_deep2_1)
ii_primus_deep2_1 = ii_primus_deep2_1[jj]
jj = get_unique_sample(ii_primus_deep2_2, frac_primus_deep2_2)
ii_primus_deep2_2 = ii_primus_deep2_2[jj]
jj = get_unique_sample(ii_primus_xmm, frac_primus_xmm)
ii_primus_xmm = ii_primus_xmm[jj]
jj = get_unique_sample(ii_dr5, frac_dr5)
ii_dr5 = ii_dr5[jj]
jj = get_unique_sample(ii_lrg, frac_lrg)
ii_lrg = ii_lrg[jj]
jj = get_unique_sample(ii_2slaq, frac_2slaq)
ii_2slaq = ii_2slaq[jj]

ii_training = [ii_zcosmos, ii_cnoc2, ii_deep2, ii_vvds, ii_primus_cosmos, ii_primus_deep2_1, ii_primus_deep2_2, ii_primus_xmm, ii_dr5, ii_lrg, ii_2slaq]
ii_training = ii[ii_training]


return, ii_training

end
