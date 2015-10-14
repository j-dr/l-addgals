pro create_catalog, name, gals, halos, outname=outname, fixa=fixa, $
                    nocolor=nocolor, declimit=declimit, des=des, $
                    vista=vista, deep=deep,johnson=johnson, ubv_deep=ubv_deep,$
                    alhambra=alhambra, swap = swap, add_bcgs=add_bcgs, dim=dim,$
                    flamex=flamex, cfhtls=cfhtls, hsc=hsc, lsst=lsst, $
		    euclid=euclid, irac=irac, wise=wise, wfirst=wfirst, $
                    twomass=twomass, galz=galz, cut=cut ;noswap=noswap

IF NOT keyword_set(des) then des=0
IF NOT keyword_set(vista) then vista=0
IF NOT keyword_set(alhambra) then alhambra=0
IF NOT keyword_set(deep) then deep=0
IF NOT keyword_set(johnson) then johnson=0
IF NOT keyword_set(ubv_deep) then ubv_deep=0
IF NOT keyword_set(flamex) then flamex=0
IF NOT keyword_set(cfhtls) then cfhtls=0
IF (keyword_set(dim)) then dim = 1
IF NOT (keyword_set(dim)) then dim = 0
IF NOT (keyword_set(cut)) then cut=0

;if not keyword_set(name) then path='' else 
;path='/data/risa/sdss/galcats/'+name+'/'
;path='/nfs/slac/g/ki/ki01/mbusha/data/sdss/galcats/'+name+'/'
path='../hv_output/'

;;;path='/nfs/slac/g/ki/ki01/mbusha/projects/addgals/sdss/galcats/'+name+'/'
if keyword_set(outname) then name = outname
;newfname=name+'f.fit' 
;;;added by mbusha to get colors working
;;;mbusha then commented all plots to let this work on orange w/ no x
;;;simpctable

;create the color catalog, if it's not already created
print, 'creating colors'
  if not keyword_set(nocolor) then assign_colors, path=path, matches=matches,$ 
    declimit=declimit,des=des, vista=vista, deep=deep, johnson=johnson, $
    ubv_deep=ubv_deep, dim=dim, alhambra=alhambra, flamex=flamex, cfhtls=cfhtls, $
    hsc=hsc, lsst=lsst, euclid=euclid, wfirst=wfirst, twomass=twomass, irac=irac, $
    wise=wise, cut=cut

print, 'reading with new file reader'
;risa_file_reader_5, gals, pathname=path
new_file_reader, gals, pathname=path, des=des, vista=vista, deep=deep,$
                 johnson=johnson, ubv_deep=ubv_deep,dim=dim, alhambra=alhambra,$
                 flamex=flamex, cfhtls=cfhtls, hsc=hsc, lsst=lsst, $
	         euclid=euclid, wfirst=wfirst, twomass=twomass, irac=irac, $
		 wise=wise, galz=galz, cut=cut
if (cut eq 0) then begin
   rdfloat,path+'gal_ginfo25.dat', amag0, amag1, amag2, amag3, amag4
endif else begin
   rdfloat,path+'gal_ginfo25.dat.cut', amag0, amag1, amag2, amag3, amag4
endelse
ss1 = create_struct('amag', fltarr(5))
amag25 = replicate(ss1, N_ELEMENTS(gals)) 
amag25.amag(0) = amag0
amag25.amag(1) = amag1
amag25.amag(2) = amag2
amag25.amag(3) = amag3
amag25.amag(4) = amag4

;ng = N_ELEMENTS(gals)
;file = '../../analysis/group_info/dataoutput/halo_fb.200.all'
;rdfloat, file, tid, tm200, tr200, trvir, x,y,z,vx,vy,vz
;for i = 0L, ng - 1 do begin
;  if (gals(i).haloid eq 527214) then continue
;  gals(i).rhalo = sqrt((gals(i).px - x(gals(i).haloid))^2 + (gals(i).py - y(gals(i).haloid))^2 + (gals(i).pz - z(gals(i).haloid))^2)
;endfor

;oplot, matches.s5, gals.dg
;Print, 'done'
;help, gals, /str
range, gals.amag(2)
;edge, gals
;this shouldn't be a problem...
;ii = where(gals.amag(2) lt -15)
;gals = gals[ii]

;ii = where(gals.dec lt 30)
;gals = gals[ii]

;plot, gals.z, gals.amag(2), psym=3, /ysty, yrange=[-24,-19]
;if not keyword_set(unswap) then begin
;unswgals = gals
if keyword_set(fixa) then fix_angles, gals
;mwrfits, gals, newfname, /create
;end

;pro swap, gals, name=name 

print, 'not swapping'
;add_brightest, gals, swgals
;******************************************************************
; mbusha changed this to remove central galaxy when adding bcg
; and re-do the galaxy ids
;******************************************************************
;IF keyword_set(add_bcgs) then begin
;IF NOT keyword_set(dim) then begin
;  central_ind = where(gals.central eq 1)
;  create_bcgs, name, mybcgs
;  centrals_to_remove = lonarr(N_ELEMENTS(mybcgs))
;  N_Removed = 0L
;  for i = 0L, N_ELEMENTS(mybcgs)-1 do begin
;    remove_ind = where(gals(central_ind).haloid eq mybcgs(i).haloid)
;    if (remove_ind(0) eq -1) then continue
;    centrals_to_remove(N_Removed) = central_ind(remove_ind(0))
;    N_Removed++
;  endfor
;  GalaxiesToKeep = lindgen(N_ELEMENTS(gals))
;  for i = 0L, N_Removed - 1 do begin
;    GalaxiesToKeep = where(GalaxiesToKeep ne centrals_to_remove(i))
;  endfor
;  gals = gals[GalaxiesToKeep]
;  gals = [mybcgs, gals]
;  for i = 0L, N_ELEMENTS(gals) - 1 do begin
;    gals(i).id = i
;  endfor
;ENDIF
;ENDIF

;IF NOT keyword_set(noswap) THEN begin
IF keyword_set(swap) THEN begin
IF NOT keyword_set(dim) THEN begin
print, 'swapping'
swap_brightest,gals,swgals
gals = swgals
endif
endif



;range, gals.amag(2)
;range, swgals.amag(2)

;swgals = swgals[ii]
;cswgals = cswgals[ii]

;help, swgals
;help, unswgals
;help, cswgals

;fname=name+'_swapped.fit' 
;unfname=name+'_unswapped.fit' 


;ii = where(gals.rhalo lt gals.r200)
;gals[ii].haloid = -1

;cswfname=name+'_galaxies.fit' 

;;;mbusha put this later to trim galaxies after we have halo info
;;;cswfname=name+'_centadd_galaxies.fit' 
;;;mwrfits, gals, cswfname, /create

;;;nname = name+'_halos.dat'
nname = name+'halos'
IF keyword_set(dim) then nname += '_dim'
nname += '.dat'
IF keyword_set(cut) then nname+='.cut'
oname = nname
;nname = '_'+nname
if not keyword_set(fixa) then begin
;;;command = '\cp '+name+'/'+'halos.dat '+nname
command = '\cp '+name+'/'+oname+' '+nname
endif else begin
;;;command = "awk '{print $1, $2, $3, $5, $4, $6}' < "+name+'/'+'halos.dat > '+nname 
command = "awk '{print $1, $2, $3, $5, $4, $6}' < "+name+'/'+oname+' > '+nname 
endelse
print, command
spawn, command

;;;hfname=name+'_halos.fit' 
hfname=name+'_halos' 
IF keyword_set(dim) then hfname += '_dim'
hfname += '.fit'
IF keyword_set(cut) then hfname+='.cut'
IF NOT keyword_set(noswap) THEN BEGIN
masstolight, gals, halos
edgehalos, halos
;mbusha commented out b/c was screwing with Millennium and we're not
;using it anyway????
;sim_ridgeline_ngals, gals, halos, nhalos
endif

;  mbusha changed to trim output *************************
;;;mwrfits, nhalos, hfname, /create
;;;cswfname=name+'_galaxies.fit' 
cswfname=name+'_galaxies' 
IF keyword_set(dim) then cswfname += '_dim'
cswfname += '.fit'
IF keyword_set(cut) then cswfname+='.cut'
gals_out = trim_gals(gals, halos)
mwrfits, gals_out, cswfname, /create

mwrfits, amag25, name+'_sdss25.fit', /create
IF (keyword_set(galz)) then begin
   galzname = name+'_galz'
   IF keyword_set(dim) then galzname += '_dim'
   galzname += '.fit'
   IF keyword_set(cut) then galzname+='.cut'
   print, 'creating struct'
   gals_galz1 = create_struct('omag',fltarr(5), 'amag',fltarr(5), 'galz', 0.0, 'pz', 0.0, 'zbin', 0L, 'central', 0L)
   print, 'replicating struct'
   gals_galz = replicate(gals_galz1,N_ELEMENTS(gals))
   print, 'setting values'
   gals_galz.omag = gals.omag
   gals_galz.amag = gals.amag
   gals_galz.galz = gals.gz
   gals_galz.pz = gals.z
   gals_galz.zbin = gals.zbin
   gals_galz.central = gals.gzcentral
   print, 'writing'
   mwrfits, gals_galz, galzname, /create
ENDIF
IF (keyword_set(des)) then begin
  desname=name+'_des'
  IF keyword_set(dim) then desname += '_dim'
  desname += '.fit'
  IF keyword_set(cut) then desname+='.cut'
  gals_des1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_des = replicate(gals_des1,N_ELEMENTS(gals))
  gals_des.omag = gals.odesmag
  gals_des.amag = gals.adesmag
  mwrfits, gals_des, desname, /create
ENDIF
IF keyword_set(vista) then begin
  vistaname=name+'_vista'
  IF keyword_set(dim) then vistaname += '_dim'
  vistaname += '.fit'
  IF keyword_set(cut) then vistaname+='.cut'
  gals_vista1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_vista = replicate(gals_vista1,N_ELEMENTS(gals))
  gals_vista.omag = gals.ovistamag
  gals_vista.amag = gals.avistamag
  mwrfits, gals_vista, vistaname, /create
ENDIF
IF keyword_set(deep) then begin
  deepname=name+'_deep'
  IF keyword_set(dim) then deepname += '_dim'
  deepname += '.fit'
  IF keyword_set(cut) then deepname+='.cut'
  gals_deep1 = create_struct('omag',fltarr(3), 'amag',fltarr(3))
  gals_deep = replicate(gals_deep1,N_ELEMENTS(gals))
  gals_deep.omag = gals.odeepmag
  gals_deep.amag = gals.adeepmag
  mwrfits, gals_deep, deepname, /create
ENDIF
IF keyword_set(johnson) then begin
  johnsonname=name+'_johnson'
  IF keyword_set(dim) then johnsonname += '_dim'
  johnsonname += '.fit'
  IF keyword_set(cut) then johnsonname+='.cut'
  gals_johnson1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_johnson = replicate(gals_johnson1,N_ELEMENTS(gals))
  gals_johnson.omag = gals.ojohnsonmag
  gals_johnson.amag = gals.ajohnsonmag
  mwrfits, gals_johnson, johnsonname, /create
ENDIF
IF keyword_set(ubv_deep) then begin
  ubv_deepname=name+'_ubv_deep'
  IF keyword_set(dim) then ubv_deepname += '_dim'
  ubv_deepname += '.fit'
  IF keyword_set(cut) then uvb_deepname+='.cut'
  gals_ubv_deep1 = create_struct('omag',fltarr(3), 'amag',fltarr(3))
  gals_ubv_deep = replicate(gals_ubv_deep1,N_ELEMENTS(gals))
  gals_ubv_deep.omag = gals.oubv_deepmag
  gals_ubv_deep.amag = gals.aubv_deepmag
  mwrfits, gals_ubv_deep, ubv_deepname, /create
ENDIF
IF keyword_set(alhambra) then begin
  alhambraname=name+'_alhambra'
  IF keyword_set(dim) then alhambraname += '_dim'
  alhambraname += '.fit'
  IF keyword_set(cut) then alhambraname+='.cut'
  gals_alhambra1 = create_struct('omag',fltarr(20), 'amag',fltarr(20))
  gals_alhambra = replicate(gals_alhambra1,N_ELEMENTS(gals))
  gals_alhambra.omag = gals.oalhambramag
  gals_alhambra.amag = gals.aalhambramag
  mwrfits, gals_alhambra, alhambraname, /create
ENDIF
IF keyword_set(flamex) then begin
  flamexname=name+'_flamex'
  IF keyword_set(dim) then flamexname += '_dim'
  flamexname += '.fit'
  IF keyword_set(cut) then flamexname+='.cut'
  gals_flamex1 = create_struct('omag',fltarr(8), 'amag',fltarr(8))
  gals_flamex = replicate(gals_flamex1,N_ELEMENTS(gals))
  gals_flamex.omag = gals.oflamexmag
  gals_flamex.amag = gals.aflamexmag
  mwrfits, gals_flamex, flamexname, /create
ENDIF
IF keyword_set(cfhtls) then begin
  cfhtlsname=name+'_cfhtls'
  IF keyword_set(dim) then cfhtlsname += '_dim'
  cfhtlsname += '.fit'
  IF keyword_set(cut) then cfhtlsname+='.cut'
  gals_cfhtls1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_cfhtls = replicate(gals_cfhtls1,N_ELEMENTS(gals))
  gals_cfhtls.omag = gals.ocfhtlsmag
  gals_cfhtls.amag = gals.acfhtlsmag
  mwrfits, gals_cfhtls, cfhtlsname, /create
ENDIF
IF keyword_set(hsc) then begin
  hscname=name+'_hsc'
  IF keyword_set(dim) then hscname += '_dim'
  hscname += '.fit'
  IF keyword_set(cut) then hscname+='.cut'
  gals_hsc1 = create_struct('omag',fltarr(5), 'amag',fltarr(5))
  gals_hsc = replicate(gals_hsc1,N_ELEMENTS(gals))
  gals_hsc.omag = gals.ohscmag
  gals_hsc.amag = gals.ahscmag
  mwrfits, gals_hsc, hscname, /create
ENDIF
IF keyword_set(lsst) then begin
  lsstname=name+'_lsst'
  IF keyword_set(dim) then lsstname += '_dim'
  lsstname += '.fit'
  IF keyword_set(cut) then lsstname+='.cut'
  gals_lsst1 = create_struct('omag',fltarr(8), 'amag',fltarr(8))
  gals_lsst = replicate(gals_lsst1,N_ELEMENTS(gals))
  gals_lsst.omag = gals.olsstmag
  gals_lsst.amag = gals.alsstmag
  mwrfits, gals_lsst, lsstname, /create
ENDIF
IF keyword_set(euclid) then begin
  euclidname=name+'_euclid'
  IF keyword_set(dim) then euclidname += '_dim'
  euclidname += '.fit'
  IF keyword_set(cut) then euclidnname+='.cut'
  gals_euclid1 = create_struct('omag',fltarr(4), 'amag',fltarr(4))
  gals_euclid = replicate(gals_euclid1,N_ELEMENTS(gals))
  gals_euclid.omag = gals.oeuclidmag
  gals_euclid.amag = gals.aeuclidmag
  mwrfits, gals_euclid, euclidname, /create
ENDIF
IF keyword_set(irac) then begin
  iracname=name+'_irac'
  IF keyword_set(dim) then iracname += '_dim'
  iracname += '.fit'
  IF keyword_set(cut) then iracname+='.cut'
  gals_irac1 = create_struct('omag',fltarr(4), 'amag',fltarr(4))
  gals_irac = replicate(gals_irac1,N_ELEMENTS(gals))
  gals_irac.omag = gals.oiracmag
  gals_irac.amag = gals.airacmag
  mwrfits, gals_irac, iracname, /create
ENDIF
IF keyword_set(wise) then begin
  wisename=name+'_wise'
  IF keyword_set(dim) then wisename += '_dim'
  wisename += '.fit'
  IF keyword_set(cut) then wisename+='.cut'
  gals_wise1 = create_struct('omag',fltarr(4), 'amag',fltarr(4))
  gals_wise = replicate(gals_wise1,N_ELEMENTS(gals))
  gals_wise.omag = gals.owisemag
  gals_wise.amag = gals.awisemag
  mwrfits, gals_wise, wisename, /create
ENDIF
IF keyword_set(wfirst) then begin
  wfirstname=name+'_wfirst'
  IF keyword_set(dim) then wfirstname += '_dim'
  wfirstname += '.fit'
  IF keyword_set(cut) then wfirstname+='.cut'
  gals_wfirst1 = create_struct('omag',fltarr(4), 'amag',fltarr(4))
  gals_wfirst = replicate(gals_wfirst1,N_ELEMENTS(gals))
  gals_wfirst.omag = gals.owfirstmag
  gals_wfirst.amag = gals.awfirstmag
  mwrfits, gals_wfirst, wfirstname, /create
ENDIF
IF keyword_set(twomass) then begin
  twomassname=name+'_twomass'
  IF keyword_set(dim) then twomassname += '_dim'
  twomassname += '.fit'
  IF keyword_set(cut) then twomassname+='.cut'
  gals_twomass1 = create_struct('omag',fltarr(3), 'amag',fltarr(3))
  gals_twomass = replicate(gals_twomass1,N_ELEMENTS(gals))
  gals_twomass.omag = gals.otwomassmag
  gals_twomass.amag = gals.atwomassmag
  mwrfits, gals_twomass, twomassname, /create
ENDIF


halos(*).edge = -99
mwrfits, halos, hfname, /create
IF NOT keyword_set(noswap) THEN BEGIN
print, 'not brightening BCGs'
;brighten_bcgs, name
endif
;;;IF keyword_set(des) THEN des_cleanup, name=name


print, "Finished."

end

PRO des_cleanup, name=name;, g, g2
hfname = name+"_galaxies.fit"
IF keyword_set(cut) then hfname+='.cut'
g = mrdfits(hfname, 1)
print, 'starting with ', n_elements(g)
minra = floor(medscat(g.ra, 0.01))
maxra = ceil(medscat(g.ra, 0.99))
mindec = floor(medscat(g.dec, 0.01))
maxdec = ceil(medscat(g.dec, 0.99))
print, minra, maxra, mindec, maxdec
ii = where(g.ra GT minra AND g.ra LT maxra AND g.dec GT mindec AND g.dec LT maxdec)
print, n_elements(g)-n_elements(ii), 'bad positions'
ii = where(g.ra GT minra AND g.ra LT maxra AND g.dec GT mindec AND g.dec LT maxdec AND g.omag[2] LT 25.3)
print, 1.0*n_elements(ii)/n_elements(g), 'percent kept after mag cut'
g2 = g[ii]
print, 'finishing with ', n_elements(g2)
;mwrfits, g2, hfname, /create
end


pro adddens, g, d1, d8

 s1=create_struct('id',0L,'omag',fltarr(5),$
		    'amag',fltarr(5),'z',0.0,'ra',0.0,'dec',0.0,$
                    'px',0.0,'py',0.0,'pz',0.0,'vx',0.0,'vy',0.0,'vz',0.0,$
                    'edge', 0, 'siglos', 0.0, $
                    'halopx', 0.0, 'halopy', 0.0, 'halopz', 0.0, $
                    'halovx', 0.0, 'halovy', 0.0, 'halovz', 0.0, $
                    'haloid',0L,'m200',0.0,'ngals',0,'r200',0.0,$
		    'rhalo',0.0,'halora',0.0,'halodec',0.0,'haloz',0.0, 'dl4', 0.0, 'd1',0.0, 'd8',0.0)
   str=replicate(s1,nobj)
   str.id = lindgen(nobj)

   input=''
   for i=0L,nobj-1 do begin

   end

return
end


;pro fix_angles, g
;
;dec  = g.ra
;ra = g.dec
;ii = where(ra lt 0)
;if (ii[0] ne -1) then ra[ii] = ra[ii]+360.
;g.ra = ra
;g.dec = dec
;
;
;dec  = g.halora
;ra = g.halodec
;if (ii[0] ne -1) then ra[ii] = ra[ii]+360.
;g.halora = ra
;g.halodec = dec
;end  



