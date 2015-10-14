PRO sim_ridgeline_ngals,simin,haloin,redhalo, redgals1, redgals2

;;5/30/06
  ;;BPK: This counts "ngals_r200" for a list of halos, and also ngals_1mpc 
  ;;     Right now, the Ngals_1mpc is not right for things with r200 < 1 
  ;;     because I use the galaxy catalog identificaiton of halo members
  ;;     which extends to r200. So for small halos with r200 < 1, the Ngals
  ;;     1 Mpc will be exactly the same as Ngals_r200.
  IF n_params() EQ 0 THEN BEGIN 
      print,'syntax: sim_ridgeline_ngals,sim,halosin,redhalos'
      return
  ENDIF 
  
;kc=mrdfits('/sdss11/data0/bens_maxBCG/jims_kcorrections_080904.fit',1)
;kc=mrdfits('~/Sites/data/sdss/jims_kcorrections_080904.fit',1)
kc=mrdfits('~/projects/addgals/idl/jims_kcorrections_080904.fit',1)
    
halos=simin[rem_dup(simin.haloid)]
halos=halos[sort(halos.z)]
nhalos=n_elements(halos)

print,ntostr(n_elements(halos))+' halos have members in the galaxy catalog'

simin=simin[sort(simin.haloz)]

;;sort the halos by z for the spline thing. this is useful for getting 0.4L*
;;for continuous z from the file of 0.4 L*'s I have.
himag=spline(kc.z,kc.lim_i,simin.haloz)

best_gr=simin.haloz*3.149+0.625 ;;these are the definitions I use in cluster finding.
best_ri=simin.haloz*0.769+0.3469

;;By construction in the mocks, the brightest thing in the halo really 
;;is the BCG. But one or
;;two could have weird spectra, such that in i-band it's not the brightest,
;;even though it is in r-band (unlikely, but possible I guess). 
;;So beware of this when you start looking at
;;individual halos. Here, I'm assuming the BCG is the brightest i-band object in the
;;halo, so I don't explicity require the galaxies to have i-band mag dimmer
;;than the BCG.

sr=where(simin.omag[1]-simin.omag[2] LT best_gr+0.1 AND $
         simin.omag[1]-simin.omag[2] GT best_gr-0.1 AND $
         simin.omag[2]-simin.omag[3] LT best_ri+0.12 AND $
         simin.omag[2]-simin.omag[3] GT best_ri-0.12 AND $
         simin.rhalo LE simin.r200 AND $
         simin.omag[3] LE himag)

cent = where(simin.omag[1]-simin.omag[2] LT best_gr+0.1 AND $
         simin.omag[1]-simin.omag[2] GT best_gr-0.1 AND $
         simin.omag[2]-simin.omag[3] LT best_ri+0.12 AND $
         simin.omag[2]-simin.omag[3] GT best_ri-0.12 AND $
         simin.rhalo LE simin.r200 AND simin.central)

sr1=where(simin.omag[1]-simin.omag[2] LT best_gr+0.1 AND $
          simin.omag[1]-simin.omag[2] GT best_gr-0.1 AND $
          simin.omag[2]-simin.omag[3] LT best_ri+0.12 AND $
          simin.omag[2]-simin.omag[3] GT best_ri-0.12 AND $
          simin.rhalo LE 1 AND $
          simin.omag[3] LE himag)


sr2=where(simin.omag[1]-simin.omag[2] LT best_gr+0.1 AND $
          simin.omag[1]-simin.omag[2] GT best_gr-0.1 AND $
          simin.omag[2]-simin.omag[3] LT best_ri+0.12 AND $
          simin.omag[2]-simin.omag[3] GT best_ri-0.12 AND $
          simin.rhalo LE 1 AND $
          simin.amag[3] LE -20.25)

sr3=where(simin.omag[1]-simin.omag[2] LT best_gr+0.1 AND $
          simin.omag[1]-simin.omag[2] GT best_gr-0.1 AND $
          simin.omag[2]-simin.omag[3] LT best_ri+0.12 AND $
          simin.omag[2]-simin.omag[3] GT best_ri-0.12 AND $
          simin.rhalo LE simin.r200 AND $
          simin.amag[3] LE -20.25-1.62*(simin.z-0.1))


sr4=where(simin.omag[1]-simin.omag[2] LT best_gr+0.1 AND $
          simin.omag[1]-simin.omag[2] GT best_gr-0.1 AND $
          simin.omag[2]-simin.omag[3] LT best_ri+0.12 AND $
          simin.omag[2]-simin.omag[3] GT best_ri-0.12 AND $
          simin.rhalo LE simin.r200)

nhalos=n_elements(haloin)

;;now use reverse indices for speed.
revh=histogram(simin[sr].haloid,rev=rev,min=0)
revh=histogram(simin[sr1].haloid,rev=rev1,min=0)
revh=histogram(simin[sr2].haloid,rev=rev2,min=0)
revh=histogram(simin[sr3].haloid,rev=rev3,min=0)
revh=histogram(simin[sr4].haloid,rev=rev4,min=0)

;add_tags,haloin,['N_ridgeline_1Mpc','N_ridgeline_r200'],$
;                ['0','0'],redhalo

add_tags,haloin,['N_ridgeline_1Mpc','N_ridgeline_r200', 'N_ridgeline0', 'N_ridgeline1', 'N_ridgeline2'],$
               ['0','0', '0', '0', '0'],redhalo

;add_tags,haloin,['N_ridgeline0', 'N_ridgeline1', 'N_ridgeline2'],$
;                [ '0', '0', '0'],redhalo


;redhalo = haloin
add_tags, simin[sr], ['N_ridgeline_r200'],$
          ['0'],redgals1

add_tags, simin[sr1], ['N_ridgeline_1mpc'],$
          ['0'],redgals2

FOR i=0L,nhalos-1 DO BEGIN 

    id = haloin[i].haloid

    ;;first the r200 cut
    IF rev[id] LT rev[id+1] THEN BEGIN 
        i1=rev[rev[id]:rev[id+1]-1]
        redhalo[i].n_ridgeline_r200=n_elements(i1)
        redgals1[i1].n_ridgeline_r200 = redhalo[i].n_ridgeline_r200 
    ENDIF 

    ;;second, the 1 Mpc cut.
    IF rev1[id] LT rev1[id+1] THEN BEGIN 
        i1=rev1[rev1[id]:rev1[id+1]-1]
        redhalo[i].n_ridgeline_1mpc=n_elements(i1)
        redgals2[i1].n_ridgeline_1mpc = redhalo[i].n_ridgeline_1mpc    
    ENDIF 

    
    IF rev2[id] LT rev2[id+1] THEN BEGIN 
        i1=rev2[rev2[id]:rev2[id+1]-1]
        redhalo[i].n_ridgeline0=n_elements(i1)
    ENDIF 

    IF rev3[id] LT rev3[id+1] THEN BEGIN 
        i1=rev3[rev3[id]:rev3[id+1]-1]
        redhalo[i].n_ridgeline1=n_elements(i1)
    ENDIF 

    IF rev4[id] LT rev4[id+1] THEN BEGIN 
        i1=rev4[rev4[id]:rev4[id+1]-1]
        redhalo[i].n_ridgeline2=n_elements(i1)
    ENDIF 

ENDFOR 
END 
