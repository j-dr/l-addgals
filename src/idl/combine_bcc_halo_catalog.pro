pro combine_bcc_halo_catalog, f1, f2, f3, outfile
;pro combine_bcc_halo_catalog, path

;outfile = path+'/Aardvark_halos.fit'
;f1 = path+'/PO_Aardvark_1050_halos.fit'
;f2 = path+'/PO_Aardvark_2600_halos.fit'
;f3 = path+'/PO_Aardvark_4000_halos.fit'

z0 = 0.0
z1 = 0.34
z2 = 0.9
z3 = 2.0

h1 = mrdfits(f1, 1)
h2 = mrdfits(f2, 1)
h3 = mrdfits(f3, 1)

ii1 = where(h1.z ge z0 and h1.z lt z1)
ii2 = where(h2.z ge z1 and h2.z lt z2)
ii3 = where(h3.z ge z2 and h3.z lt z3)

h = [h1[ii1], h2[ii2], h3[ii3]]
h.ra = atan2(h.halopy, h.halopx, /deg)

mwrfits, h, outfile, /create

end

