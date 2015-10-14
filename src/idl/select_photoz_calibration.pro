pro select_photoz_calibration,gal,ind,npick=npick,area=area, ra_range=ra_range

   if not keyword_set(npick) then npick=20000.
   if not keyword_set(area) then area=573.
   if not keyword_set(ra_range) then ra_range=[0.,10.]
   ra_min = 2.5+ra_range(0)
   ra_width = (ra_range(1)-ra_range(0))-5
   ntotal=n_elements(gal)


   sdss=where(gal.omag(1) le 17.8)

   deep=where(gal.omag(3) le 23)
   ndeep=n_elements(deep)
   npick_deep=54000.
;   area_search = (npick_deep/ndeep)*225.
   area_search = (npick_deep/ndeep)*area
   radius = sqrt(area_search/!pi)
;   random_center_ra = 2.5+randomu(seed)*10.
;   random_center_dec = 2.5+randomu(seed)*10.
   random_center_ra = ra_min+randomu(seed)*ra_width
   random_center_dec = ra_min+randomu(seed)*ra_width
   close_match_radec,random_center_ra,random_center_dec,$
     gal(deep).ra,gal(deep).dec,i1,i2,radius,100000
   deep=deep(i2)

   vvdss=where(gal.omag(3) le 22.5)
   nvvdss=n_elements(vvdss)
   npick_vvdss=26000.
;   area_search = (26000./nvvdss)*225.
   area_search = (26000./nvvdss)*area
   radius = sqrt(area_search/!pi)
;   random_center_ra = 2.5+randomu(seed)*10.
;   random_center_dec = 2.5+randomu(seed)*10.
   random_center_ra = ra_min+randomu(seed)*ra_width
   random_center_dec = ra_min+randomu(seed)*ra_width
   close_match_radec,random_center_ra,random_center_dec,$
     gal(vvdss).ra,gal(vvdss).dec,i1,i2,radius,100000
   vvdss=vvdss(i2)

   vvdsd=where(gal.omag(3) le 24)
   nvvdsd=n_elements(vvdsd)
   npick_vvdsd=60000.
;   area_search = (60000./nvvdsd)*225.
   area_search = (60000./nvvdsd)*area
   radius = sqrt(area_search/!pi)
;   random_center_ra = 2.5+randomu(seed)*10.
;   random_center_dec = 2.5+randomu(seed)*10.
   random_center_ra = ra_min+randomu(seed)*ra_width
   random_center_dec = ra_min+randomu(seed)*ra_width
   close_match_radec,random_center_ra,random_center_dec,$
     gal(vvdsd).ra,gal(vvdsd).dec,i1,i2,radius,100000
   vvdsd=vvdsd(i2)

   ind = [sdss,deep,vvdss, vvdsd]
 
   return
   end
