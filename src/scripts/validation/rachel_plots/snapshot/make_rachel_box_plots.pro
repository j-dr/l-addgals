pro make_rachel_box_plots, dirlist, labels, outdir, typelist=typelist, wplist=wplist

;general, flexible plotting set
;dir = "/nfs/slac/g/ki/ki11/des/mbusha/catalogs/Consuelo/addgals_test/plot_info/"
if not KEYWORD_SET(typelist) then typelist = "consuelo"
if not KEYWORD_SET(wplist) then wplist = "sham"

spawn,"mkdir -p "+outdir
dir = ''
;labels = labels[0:4]
;dirlist = dirlist[0:4]
;typelist = typelist[0:4]
;wplist = wplist[0:4]

nobar=1

;LF
infile = "~rmredd/ki10/make_sham/L80G/smf_mcmc2/data/LF_new.dat"
;infile = "/nfs/slac/g/ki/ki11/des/mbusha/catalogs/Consuelo/addgals_test/v05/LF.dat"
plot_addgals_lf,outdir,infile,files=dir+dirlist+"/"+typelist+"_lf.idlsav",$
  labels=labels,nobar=nobar

;HOD
plot_gen_hod,outdir+"hod_addgals_testing.ps",dir+dirlist+$
  "/"+typelist+"_hod.idlsav",labels=labels,nobar=nobar

;Mhost(Mr)
;plot_mr_mhost,outdir,dir+dirlist+"/"+typelist+"_mr_mhost_reln.idlsav",$
;  labels,nobar=nobar
;plot_mr_mhost,outdir+'5e12_',dir+dirlist+"/5e12_consuelo_mr_mhost_reln.idlsav",$
;  labels

;CLF
plot_gen_clfs,outdir,dir+dirlist+"/"+typelist+"_clf.idlsav",$
  labels=labels,type=0,nobar=nobar
plot_gen_clfs,outdir,dir+dirlist+"/"+typelist+"_clf.idlsav",$
  labels=labels,type=1,nobar=nobar
plot_gen_clfs,outdir,dir+dirlist+"/"+typelist+"_clf.idlsav",$
  labels=labels,type=2,nobar=nobar
plot_clf_sum,outdir,dir+dirlist+"/"+typelist+"_clf.idlsav",labels=labels

;rpr
plot_gen_rpr,outdir+"rpr_19_addgals_testing_cons.ps",$
  dir+dirlist+"/"+typelist+"_rpr_proj_19.idlsav",labels=labels,/full,/nosc,nobar=nobar
plot_gen_rpr,outdir+"rpr_20_addgals_testing_cons.ps",$
  dir+dirlist+"/"+typelist+"_rpr_proj_20.idlsav",labels=labels,/full,/nosc,nobar=nobar
plot_gen_rpr,outdir+"rpr_21_addgals_testing_cons.ps",$
  dir+dirlist+"/"+typelist+"_rpr_proj_21.idlsav",labels=labels,/full,/nosc,nobar=nobar

;wp
set_plot,'x'
multiplot,/reset
plot_gen_wp_set,outdir,filenames=dir+dirlist+"/"+typelist+"_"+wplist+"_wp_mr_",$
  labels=labels,nobar=nobar

;P(B!=C)

;fsat
;plot_add_fsat,outdir,dir+dirlist+"/"+typelist+"_fsat_m12.7.idlsav",labels

end
