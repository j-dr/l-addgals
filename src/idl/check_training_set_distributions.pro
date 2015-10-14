skip_read = 1

if (skip_read eq 0) then begin
  ;;;out input training set and SDSS magnitudes
  in='/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.5d/photoz_DR8/PO_Aardvark_DR8_training_set.fit'
  in2='/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.5d/photoz_DR8/PO_Aardvark_DR8_training_set_sdss_mag.fit'
  t = mrdfits(in, 1)
  tmag = mrdfits(in2, 1)
  
  ;;;out input photometric set
  pin1='/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.5d/truth_hod/Aardvark_v0.5d_truth_des_masked.97.fit'
  pin12='/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.5d/DR8/Aardvark_v0.5d_sdss_mag.97.fit'
  g1 = mrdfits(pin1,1)
  gmag1 = mrdfits(pin12,1)
  pin2='/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.5d/truth_hod/Aardvark_v0.5d_truth_des_masked.86.fit'
  pin22='/nfs/slac/g/ki/ki18/des/mbusha/catalogs/Aardvark/Catalog_v0.5d/DR8/Aardvark_v0.5d_sdss_mag.86.fit'
  g2 = mrdfits(pin2,1)
  gmag2 = mrdfits(pin22,1)
  g = [g1, g2]
  gmag = [gmag1, gmag2]
endif

;;;select objects that go to 21.8 with no color selection
ii_zcosmos =        where(t.ra gt 140 and t.ra lt 141  and t.dec gt 41   and t.dec lt 43.3, n_zcosmos)
ii_cnoc2 =          where(t.ra gt 61 and t.ra lt 62.145    and t.dec gt 29   and t.dec lt 29.4, n_conc2)
ii_deep2 =          where(t.ra gt 129 and t.ra lt 130.14    and t.dec gt 29   and t.dec lt 29.4, n_deep2)
ii_vvds =           where(t.ra gt 10 and t.ra lt 12.2    and t.dec gt 23   and t.dec lt 25, n_vvds)
ii_primus_cosmos =  where(t.ra gt 62 and t.ra lt 63 and t.dec gt 18   and t.dec lt 19.16, n_primus_c)
ii_primus_deep2_1 = where(t.ra gt 33 and t.ra lt 34 and t.dec gt 14   and t.dec lt 14.72, n_primus_d1)
ii_primus_deep2_2 = where(t.ra gt 151 and t.ra lt 152    and t.dec gt 9 and t.dec lt 9.71, n_primus_d2)
ii_primus_xmm =     where(t.ra gt 157 and t.ra lt 158.71    and t.dec gt 4   and t.dec lt 6, n_primus_xmm)

iitrain = [ii_zcosmos, ii_cnoc2, ii_deep2, ii_vvds, ii_primus_cosmos, ii_primus_deep2_1, ii_primus_deep2_2, ii_primus_xmm]

;;;select the objects in the trianing set that come from DR5, they're close to the end
target_2slaq = 8633
target_lrg = 20000L
target_dr5 = 100000L
idr5 = lindgen(target_dr5)+N_ELEMENTS(t) - target_dr5-target_lrg-target_2slaq

;;;the dr5 objects are cut to r-band 17.8.  Make the same galaxy cut
ig = where(gmag.omag[2] lt 17.8)

;;;compare the magnitude distributions in the g, r, i bands
begplot, '~/temp/training_photo_mags.ps', /color, xsize = 7, ysize = 6
mbin = 0.1
pplothist, tmag[idr5].omag[2], bin = mbin, /prob, /ylog, xrange=[12,19], $
	yrange = [1e-3,10], xtitle = 'mag', ytitle = 'P(mag)'
pplothist, tmag[idr5].omag[2], bin = mbin, /prob, /over, color = !red
pplothist, gmag[ig].omag[2], bin = mbin, /over, color=!red, /prob,linestyle=2
pplothist, tmag[idr5].omag[1], bin = mbin, /prob, /over, color = !green
pplothist, gmag[ig].omag[1], bin = mbin, /over, color=!green, /prob,linestyle=2
pplothist, tmag[idr5].omag[3], bin = mbin, /prob, /over, color = !blue
pplothist, gmag[ig].omag[3], bin = mbin, /over, color=!blue, /prob,linestyle=2
legend, ['m!Dg', 'm!Dr', 'm!Di'], linestyle = 0, color = [!green, !red, !blue]
endplot


;;;compare the redshift distributions
begplot, '~/temp/training_photo_z.ps', /color, xsize = 7, ysize = 6
zbin = 0.03
pplothist, t[idr5].z, bin = zbin, /prob, xtitle = 'z', ytitle = 'P(z)'
pplothist, g[ig].z, bin = zbin, /prob, /over, color = !red
legend, ['DR5 Training', 'Photometric'], linestyle=0, $
	color = [!P.COLOR, !red], /right
endplot


;;;compare to the non-dr5 stuff
gp = cut_to_dr8_training_limits(gmag, ii=iicut)
ig = where(gmag.omag[2] lt 21.8)

mbin = 0.2
pplothist, tmag[iitrain].omag[2], bin = mbin, /prob, /ylog, xrange=[13,22], $
        yrange = [1e-3,10], xtitle = 'mag', ytitle = 'P(mag)'
pplothist, tmag[iitrain].omag[2], bin = mbin, /prob, /over, color = !red
pplothist, gmag[ig].omag[2], bin = mbin, /over, color=!red, /prob,linestyle=2
pplothist, gp.omag[2], bin = mbin, /over, color=!red, /prob,linestyle=1
pplothist, tmag[iitrain].omag[1], bin = mbin, /prob, /over, color = !green
pplothist, gmag[ig].omag[1], bin = mbin, /over, color=!green, /prob,linestyle=2
pplothist, gp.omag[1], bin = mbin, /over, color=!green, /prob,linestyle=1
pplothist, tmag[iitrain].omag[3], bin = mbin, /prob, /over, color = !blue
pplothist, gmag[ig].omag[3], bin = mbin, /over, color=!blue, /prob,linestyle=2
pplothist, gp.omag[3], bin = mbin, /over, color=!blue, /prob,linestyle=1
legend, ['m!Dg', 'm!Dr', 'm!Di'], linestyle = 0, color = [!green, !red, !blue]


begplot, '~/temp/training_photo_z22.ps', /color, xsize = 7, ysize = 6
zbin = 0.06
pplothist, t[iitrain].z, bin = zbin, /prob, xtitle = 'z', ytitle = 'P(z)'
pplothist, g[ig].z, bin = zbin, /prob, /over, color = !red
pplothist, g[iicut].z, bin = zbin, /prob, /over, color = !blue
legend, ['DR5 Training', 'Photometric with r < 21.8', 'Photometric with u,r < 24.6, 21.8'], linestyle=0, $
        color = [!P.COLOR, !red, !blue], /right
endplot

end
