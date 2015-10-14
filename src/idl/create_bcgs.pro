;PRO add_bcgs, name, orig=orig, new=new, bcgs

pro plot_this, bcgs, hh, index
;g = mrdfits('hv_v1.00ja_galaxies.fit', 1)
bcgs = mrdfits('maxbcg_uncut_cluster_catalog_n200ge5.fit', 1)
name =  'hv_v1.00ja'
hfile = name+'_halos.fit'
h = mrdfits(hfile, 1)
ii = where(h.z gt 0.05 and h.z lt 0.3)
hh = h[ii]
add_bcgs, name, bcgs, index
plot, hh.m200, bcgs[index].bcg_ilum, psym=3, /xlog, /ylog
bin_ave, hh.m200, bcgs[index].bcg_ilum, bmin=13.3, n=10  
add_tags, hh, ['cmodel_counts', 'bcg_ilum', 'imag', 'num2_imag'], ['dblarr(5)', '0.0d', '0.0d', '0.0d'], hhbcgs
hhbcgs.cmodel_counts = bcgs[index].cmodel_counts
hhbcgs.bcg_ilum = bcgs[index].bcg_ilum
hhbcgs.imag = bcgs[index].imag
hhbcgs.num2_imag = bcgs[index].num2_imag
mwrfits, hhbcgs, 'halos_bcgs.fit', /create
end

pro tmp_plot, bcgs
!p.multi=[0,2, 2]
plothist, alog10(1.d10*bcgs[where(bcgs.ngals gt 10)].bcg_ilum), bin=0.01, /norm
plot, bcgs.ngals_R200, alog10(bcgs.bcg_ilum), psym=3, /xlog  

hfile = 'hv_v1.00ja_halos.fit'
h = mrdfits(hfile, 1)
mass = h.m200
lum = 5.8d6*mass^0.29
scat_mult = 0.2*(mass/1e14)^(-0.05)

scatter = scat_mult*randomn(1, n_elements(mass))
lum = 10.^(alog10(lum)+scatter)
plothist, alog10(lum), bin=0.01, /over, color=!red, /norm
plot, h.m200, scat_mult, psym=3, /xlog, /ysty
plot, h.m200, alog10(lum), psym=3, /xlog, /ysty
bin_ave, h.m200, alog10(lum), bmin=13.5, n=10, bin=0.15
oplot, h.m200, alog10(5.8d6*h.m200^0.29)
oplot, h.m200, alog10(5.8d6*h.m200^0.29)+0.2, thick=4
oplot, h.m200, alog10(5.8d6*h.m200^0.29)-0.2, thick=4

end


PRO get_bcgs, name, bcgs
!p.multi=[0,2,2]
inbcgs = mrdfits('maxbcg_uncut_cluster_catalog_n200ge3.fit', 1)
inbcgs = inbcgs[sort(inbcgs.bcg_ilum)]
;name =  'hv_v1.00ja'

;********************************;
;   readin changed by mbusha    ;
;********************************;
;;;hfile = name+'_halos.fit'
;;;h = mrdfits(hfile, 1)
hfile = name+'_halos.dat'
h = read_halo_file(hfile)

;;;minimum galue changed by mbusha to agree with sdss
;;;ii = where(h.z gt 0.05 and h.z lt 0.3)
ii = where(h.z gt 0.056 and h.z lt 0.3)
h = h[ii]
add_the_bcg, h, inbcgs, index

;add_tags, hh, ['cmodel_counts', 'bcg_ilum', 'imag', 'num2_imag'], ['dblarr(5)', '0.0d', '0.0d', '0.0d'], hhbcgs
;hhbcgs.cmodel_counts = bcgs[index].cmodel_counts
;hhbcgs.bcg_ilum = bcgs[index].bcg_ilum
;hhbcgs.imag = bcgs[index].imag
;hhbcgs.num2_imag = bcgs[index].num2_imag
;mwrfits, hhbcgs, 'halos_bcgs.fit', /create

nbcgs = n_elements(h)

;;;changed by mbusha to make consistent with other gals definition
;s1=create_struct('id',0L,'omag',fltarr(5),$
;		    'amag',fltarr(5),'z',0.0,'ra',0.0,'dec',0.0,$
;                    'px',0.0,'py',0.0,'pz',0.0,'vx',0.0,'vy',0.0,'vz',0.0,$
;                    'nonedge', 0, 'siglos', 0.0, $
;                    'halopx', 0.0, 'halopy', 0.0, 'halopz', 0.0, $
;                    'halovx', 0.0, 'halovy', 0.0, 'halovz', 0.0, $
;                    'haloid',0L,'m200',0.0,'ngals',0,'r200',0.0,$
;
;                    'rhalo',0.0,'halora',0.0,'halodec',0.0,'haloz',0.0, 'd8', 0.0, 'central',0)
s1=create_struct('id',0L,'ecatid',0L,'omag',fltarr(5),$
                 'amag',fltarr(5),'z',0.0,'ra',0.0,'dec',0.0,$
                 'px',0.0,'py',0.0,'pz',0.0,'vx',0.0,'vy',0.0,'vz',0.0,$
                 'edge', 0, 'siglos', 0.0, $
                 'halopx', 0.0, 'halopy', 0.0, 'halopz', 0.0, $
                 'halovx', 0.0, 'halovy', 0.0, 'halovz', 0.0, $
                 'haloid',0L,'m200',0.0,'ngals',0,'r200',0.0,$
                 'rhalo',0.0,'halora',0.0,'halodec',0.0,'haloz',0.0,'central',0)

bcgs=replicate(s1,nbcgs)
for i = 0, nbcgs-1 do begin
   bcgs[i].ecatid = index[i]
   bcgs[i].omag = inbcgs[index[i]].cmodel_counts
   bcgs[i].z = h[i].z
   bcgs[i].ra = h[i].ra
   bcgs[i].dec = h[i].dec
   bcgs[i].px = h[i].halopx
   bcgs[i].py = h[i].halopy
   bcgs[i].pz = h[i].halopz
   bcgs[i].vx = h[i].halovx
   bcgs[i].vy = h[i].halovy
   bcgs[i].vz = h[i].halovz
   bcgs[i].siglos = h[i].siglos
   bcgs[i].halopx = h[i].halopx
   bcgs[i].halopy = h[i].halopy
   bcgs[i].halopz = h[i].halopz
   bcgs[i].halovx = h[i].halovx
   bcgs[i].halovy = h[i].halovy
   bcgs[i].halovz = h[i].halovz
   bcgs[i].haloid = h[i].haloid
   bcgs[i].m200 = h[i].m200
   bcgs[i].r200 = h[i].r200
   bcgs[i].rhalo = 0
   bcgs[i].ngals = h[i].ngals
   bcgs[i].haloz = h[i].z
   bcgs[i].halora = h[i].ra
   bcgs[i].halodec = h[i].dec
   bcgs[i].central = 1
endfor
end


PRO add_bcgs_tohalos, name, bcgs, index
!p.multi=[0,3,2]
bcgs = mrdfits('maxbcg_uncut_cluster_catalog_n200ge3.fit', 1)
bcgs = bcgs[sort(bcgs.bcg_ilum)]

;name =  'hv_v1.00ja'
hfile = name+'_halos.fit'
h = mrdfits(hfile, 1)
ii = where(h.z gt 0.1 and h.z lt 0.11)
hh = h[ii]
help, hh
add_the_bcg, hh, bcgs, index
;plot, hh.m200, bcgs[index].bcg_ilum, psym=3, /xlog, /ylog
;bin_ave, hh.m200, bcgs[index].bcg_ilum, bmin=13.3, n=10  

plot, hh.m200, bcgs[index].bcg_ilum, psym=3, /xlog, /ysty
bin_ave, hh.m200, bcgs[index].bcg_ilum, bmin=13.5, n=10, bin=0.15

add_tags, hh, ['cmodel_counts', 'bcg_ilum', 'imag', 'num2_imag'], ['dblarr(5)', '0.0d', '0.0d', '0.0d'], hhbcgs
hhbcgs.cmodel_counts = bcgs[index].cmodel_counts
hhbcgs.bcg_ilum = bcgs[index].bcg_ilum
hhbcgs.imag = bcgs[index].imag
hhbcgs.num2_imag = bcgs[index].num2_imag
mwrfits, hhbcgs, 'halos_bcgs.fit', /create

end

PRO add_the_bcg, halos, bcgs, index
mass = halos.m200
lum = 5.8d6*mass^0.29
range, lum
;;;plothist, alog10(lum), bin=0.01
;;;plothist, alog10(1.d10*bcgs.bcg_ilum), /over, bin=0.01, color=!blue
nnn = n_elements(mass)
scatter = 0.2*randomu(1, nnn)
lum = 10.^(alog10(lum)+scatter)
;;;plothist, alog10(1.d10*bcgs.bcg_ilum), /over, color=!red, over, bin=0.01
;mag_r = -2.5*alog10(lum)+4.76
;plot, mass, lum, /xlog, /ylog, psym=3, xrange=[5e13, 5e15]

nz = indgen(nnn)
zzz = findgen(nnn)
lll = findgen(nnn)
index = replicate(0L, nnn)
for hi = 0L, nnn-1 do begin
   ii = where(bcgs.z gt halos[hi].z-0.003 and bcgs.z lt halos[hi].z+0.003)
   binary_search, 1.d10*bcgs[ii].bcg_ilum, lum[hi], close, /round, /edgedefault
zzz[hi] = bcgs[ii[close]].z
lll[hi] = 1.e10*bcgs[ii[close]].bcg_ilum
   index[hi] = ii[close]
  ; print, halos[hi].z, min(bcgs.z), n_elements(ii), lum[hi],lll[hi] 
endfor
;;;autohist, zzz-halos.z
;;;plot, lum, lum, /xlog, /ylog
;;;oplot, lum, lll, psym=2, color=!red
;;;autohist, (lum-lll)/lum
;aa = absmag(cmodel_counts[2]
;red,omega0=0.3,omegalambda=0.7,h100=1.0 
;plot, lum[0:nnn-1],  bcgs[0:nnn-1].cmodel_counts[3], psym=3, /ysty, /xlog
;plot, lum[0:nnn-1],  bcgs[index].cmodel_counts[3], psym=3, /ysty, /xlog

;niceprint, dmodulus(h[0:nnn-1].z), -2.5*alog10(lum[0:nnn-1])+4.76, bcgs[0:nnn-1].cmodel_counts[3]
;autohist,nz
;!p.multi=0

end
 
