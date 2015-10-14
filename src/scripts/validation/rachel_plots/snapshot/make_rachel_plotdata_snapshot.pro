;+
; NAME:
;	MAKE_RACHEL_PLOTDATA_SNAPSHOT
;
; PURPOSE:
;
;	Is a wrapper routine for generating Rachel's PaperI for an 
;	ADDGALS catalog run on one of the Consuelo boxes (SHAM or
;	ADDGALS).  
;
; NOTE:
;
;	In order to run on a SHAM catalog, it must first be
;	converted to a .fits ADDGALS format using, i.e., 
;	/afs/slac.stanford.edu/u/ki/mbusha/projects/modules/idl/addgals/rockstar_sham_to_fits.pro
;
; AUTHOR:
;
;	Rachel Reddick and Michael Busha
;
; CALLING SEQUENCE:
;
;	 make_rachel_plotdata_snapshot, '/path/to/output/', '/path/to/input/' 
;
; OPTIONAL INPUTS:
;
;	gfname:		set to something other than the default galaxies.fit
;	hfname:		set to something other than the default halos.fit
;
; KEYWORD_PARAMETERS:
;
;       old_hfile:      set to true to always read the halos locations 
;                       from the v16 catalog
;
; MODIFICATION HISTORY:
;
;	turned into a procedure by mtb on 29.01.2012.
;
;-

pro make_rachel_plotdata_snapshot, outdir, dir, gfname=gfname, hfname=hfname, $
	masscut=masscut

; setup the input filenames
if not KEYWORD_SET(gfname) then gfname = "galaxies.fit"
if not KEYWORD_SET(hfname) then hfname = "halos.fit"
gfile = dir+gfname
hfile = dir+hfname

; make sure the output directory exists
spawn,'mkdir -p '+outdir
print,outdir

; read the data
gc = mrdfits(gfile,1)
hc = mrdfits(hfile,1)

;setup host mass cut
if not KEYWORD_SET(masscut) then masscut = 2e11

; sanity checking in rseponse to some previous bad catalogs
if n_elements(where(hc.halopx eq 0 or hc.halopy eq 0 or hc.halopz eq 0)) gt 10 then hpifail = 1 else hpifail = 0
clist = where(gc.central eq 1)
if hpifail then print,"HPIFAIL"
if hpifail then hc.halopx = gc[clist].px
if hpifail then hc.halopy = gc[clist].py
if hpifail then hc.halopz = gc[clist].pz

;compile routines for doing all the calculations
resolve_routine, 'make_base_hod_fast', /either, /comp
resolve_routine, 'plot_addgals_clfs', /either, /comp
resolve_routine, 'make_clf_boot_mb', /either, /comp
resolve_routine, 'get_fsat', /either, /comp
resolve_routine, 'get_wp_addgals', /either, /comp
resolve_routine, 'make_rpr_mb_hlist', /either, /comp
resolve_routine, 'addgals_testing', /either, /comp

;run for LF/wp with no cuts
addgals_testing,outdir=outdir,$
  check_lf=1,$
  run_mhostmr=0,$
  run_hod=0,$
  run_clf=0,$
  run_color_clf=0,$
  run_fsat=0,$
  run_sigbcg=0,$
  run_wpthresh=1,$
  run_wp_color=0,$
  run_rpr=0,$
  run_pbnc=0,$
  noplot=1,$                    ;toggle plotting
  gc=gc,hc=hc,$                 ;consuelo input
  gb=gb, hb=hb,$                ;bolshoi input
  g,h                           ;carmen input


;take limits and do haloid reassignment
mlist = where(gc.amag[2] lt -18)
gc = gc[mlist]
if TAG_EXIST(hc[0], 'mvir') then hc.m200 = hc.mvir
if TAG_EXIST(hc[0], 'rvir') then hc.r200 = hc.rvir*1000.
if (max(hc.r200) lt 1000.) then hc.r200 *= 1000.
mlist = where(hc.m200 ge 2e11)
hc = hc[mlist]

haloid = get_fixed_hids(hc.haloid,hc.r200/1000,hc.m200,hc.halopx,hc.halopy,hc.halopz,$
                        gc.px,gc.py,gc.pz,rhalo=rhalo,central=central)

index = lonarr(max(hc.haloid+1))
index[hc.haloid] = lindgen(n_elements(hc.haloid))
gc.haloid = haloid
gc.rhalo = rhalo
gc.central = central
list = where(haloid ge 0)
gc.m200 = 0*gc.m200
gc[list].m200 = hc[index[haloid[list]]].m200
gc.r200 = 0*gc.r200
gc[list].r200 = hc[index[haloid[list]]].r200

mlist = where(gc.m200 gt masscut)
gc = gc[mlist]

;run with given mass cut
addgals_testing,outdir=outdir,$
  check_lf=0,$
  run_mhostmr=0,$
  run_hod=1,$
  run_clf=1,$
  run_color_clf=0,$
  run_fsat=1,$
  run_sigbcg=0,$
  run_wpthresh=0,$
  run_wp_color=0,$
  run_rpr=1,$
  run_pbnc=0,$
  noplot=1,$                    ;toggle plotting
  gc=gc,hc=hc,$                 ;consuelo input
  gb=gb, hb=hb,$                ;bolshoi input
  g,h                           ;carmen input

end
