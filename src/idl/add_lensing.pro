;+
; ROUTINE:  add_lensing
; PURPOSE:  lens a galaxy catalog using a shear catalog
; USEAGE:   add_lensing, g, shear
; INPUT:    
;  g - a galaxy catalog
;  shear - a shear catalog
; OUTPUT: lensed galaxy catalog in g
; NOTES:  see the comments below to understand what the code is doing
; AUTHOR:   Michael Busha based on code written by Joerg Dietrich, early 2012
; REVISIONS:
;    Swapped sign on g1 to match ra-dec convention, not DECam
;    Fixed a bug in how galaxies were lensed when |g| > 1
;    Added parameter w to capture rotation
;     -Matthew Becker, 8-30-2012
;-
;

pro add_lensing, g, shear
  
;;check to make sure all the tags we need exist 
;;assume if tra is not there, then none of tags are there

if (tag_exist(g, 'tra') eq 0) then begin
   tagnames = ['tra', 'tdec', 'gamma1', 'gamma2', 'kappa', 'w', 'mu', 'size', 'lmag', 'epsilon']
   values = ['0.', '0.', '0.', '0.', '0.', '0.','0.', '0.', 'fltarr(5)', 'fltarr(2)']
   add_tags, g, tagnames, values
   g.tra = g.ra
   g.tdec = g.dec
endif

if (tag_exist(g, 'w') eq 0) then begin
   tagnames = ['w']
   values = ['0.']
   add_tags, g, tagnames, values, ng
   g = ng
endif


;;gamma1 == 0 marks galaxies which have zero images
g.gamma1 = 0.
for i = 0L, N_ELEMENTS(shear) - 1 do begin
  ;;get the index of the source for this image
  tind = shear[i].index

  ;;if gamma1 != 0, then this galaxy is multiply imaged
  ;;thus we add it to the end of the catalog
  if (g[tind].gamma1 ne 0.) then begin 
    g = [g, g[tind]]
    tind = N_ELEMENTS(g) - 1
  endif

  g[tind].ra = shear[i].ra
  g[tind].dec = shear[i].dec
  
  ;;extract g1,g2,kappa,w from A using formulas from Vale & White (2003)
  ;; A = |1-kappa-gamma1      -gamma2-w|
  ;;     |-gamma2+w      1-kappa+gamma1|

  g[tind].gamma1 = 0.5*(shear[i].a11 - shear[i].a00)
  g[tind].gamma2 = -0.5*(shear[i].a01 + shear[i].a10)
  g[tind].kappa = 1.0 - 0.5*(shear[i].a00 + shear[i].a11)
  g[tind].w = 0.5*(shear[i].a10 - shear[i].a01)
  
  ;;compute mu = 1/detA
  g[tind].mu = 1./(shear[i].a11*shear[i].a00 - shear[i].a01*shear[i].a10)
  
  ;;lens the size and magnitudes
  g[tind].size = g[tind].tsize*sqrt(g[tind].mu)
  for im = 0, 4 do g[tind].lmag[im] = g[tind].tmag[im] - 2.5*alog10(g[tind].mu)
  
  ;;get intrinsic shape
  epss = complex(g[tind].te[0], g[tind].te[1])
  
  ;;get reduced complex shear g = (gamma1 + i*gamma2)/(1-kappa)
  gquot = complex(g[tind].gamma1, g[tind].gamma2) / complex(1.-g[tind].kappa, 0.)
  
  ;;lens the shapes - see Bartelmann & Schneider (2001), Section 4.2
  if (abs(gquot) lt 1) then begin
    eps = (epss+gquot)/(complex(1.0,0.0)+(epss*conj(gquot)))
  endif else begin
    eps = (complex(1.0, 0.0)+(gquot*conj(epss)))/(conj(epss)+conj(gquot))
  endelse
  g[tind].epsilon[0] = real_part(eps)
  g[tind].epsilon[1] = imaginary(eps)
endfor

end
