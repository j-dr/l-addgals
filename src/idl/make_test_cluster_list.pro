
if (0) then begin
for i = 0L, 1 do begin
  file = '~/ki01/projects/addgals/PO_oct/v1.07_z/idl/cluster_comp_v1.07_fixedID'+strcompress(string(i),/remove_all)+'.fits'
  g_tmp = mrdfits(file,1)
  if (i eq 0) then g = g_tmp
  if (i gt 0) then g = [g, g_tmp]
  g_tmp = 1
endfor
endif

mismatch = 0.1
dropped = 0.1
double = 0.01
;mismatch = 0.0
;dropped = 0.0
double = 0.0

ng = N_ELEMENTS(g)
nh_max = lonarr(max(g.haloid)+1)
ran = randomu(seed,ng)
ih = where(g.central)

nbins = 20
decdel = 90./nbins
max_in_bin =  10000
ih = lonarr(nbins,max_in_bin)
n_in_bin = lonarr(nbins)
for i = 0, nbins - 1 do begin
  decmin = i*decdel
  decmax = decmin + decdel
  ka = where(g.central and g.dec ge decmin and g.dec lt decmax, count)
  ih(i,0:count-1) = ka
  n_in_bin(i) = count
endfor

gout_ent = create_struct('rank', 0L, 'galaxyID', 0L, 'haloID', 0L, 'haloID2',0L, 'omag', 0., 'ra', 0., 'dec', 0., 'photoz', 0.)
gout = replicate(gout_ent, ng)
gout.GalaxyID = g.ID
gout.haloID = g.haloid
gout(*).haloID2 = -1
gout.omag = g.omag(1)
gout.ra = g.ra
gout.dec = g.dec
gout.photoz = g.photoz
theta = !PI*g.ra/180.
phi = (90.-g.dec)*!PI/180.

print, "Assigning to halos"
for i = 0L, ng - 1 do begin
  if (ran(i) gt mismatch+dropped) then continue
  dbin = floor(g(i).dec/decdel)
  angsep = fltarr(n_in_bin(dbin))
  angsep = cos(phi(i))*cos(phi(ih(dbin,0:n_in_bin(dbin)-1))) + sin(phi(i))*sin(phi(ih(dbin,0:n_in_bin(dbin)-1)))*cos(theta(i)-theta(ih(dbin,0:n_in_bin(dbin)-1)))
  angsep = acos(angsep)
  srt = sort(abs(angsep))
  if (g(i).haloid gt -0.5) then begin
    if (ran(i) lt mismatch) then begin
      num = floor(randomu(seed)*3)
      num = 0
      gout(i).haloID = g(ih(dbin,srt(num))).haloID
      if (ran(i) lt double) then begin
        num += 2    
        gout(i).haloID2 = g(ih(dbin,srt(num))).haloID
      endif
    endif else begin
      g(i).haloID = -1
    endelse
  endif else begin
    if (ran(i) gt mismatch) then continue
    num = floor(randomu(seed)*5)
    gout(i).haloID = g(ih(dbin,srt(num))).haloID
  endelse
endfor

;;;remake the galaxy structure
ngout = long(ng)
for i = 0L, ng - 1 do begin
  if (gout(i).haloID2 gt -0.5) then ngout += 1L
endfor

gals = replicate(gout_ent, ngout)
j = 0L
for i = 0L, ng - 1 do begin
  gals(j) = gout(i)
  j++
  if (gout(i).haloid2 gt -0.5) then begin
    gals(j) = gout(i)
    gals(j).haloID = gout(j).haloID2
    j++
  endif
endfor

;;;calculate halo info and rand
h_ent = create_struct('rank', 0L, 'ra', 0., 'dec', 0., 'out', 0, 'ngals', 0L, 'omag', 0., 'photoz', 0.)
h = replicate(h_ent, max(gals.haloID)+1)
h(*).omag = 100.
gals.rank = -1
for i = 0L, ngout - 1 do begin
  ih = gals(i).haloID
  if (ih lt 0) then continue
  h(ih).out = 1
  h(ih).ngals++
  h(ih).photoz += gals(i).photoz
  if (gals(i).omag lt h(ih).omag) then begin
    h(ih).ra = gals(i).ra
    h(ih).dec = gals(i).dec
    h(ih).omag = gals(i).omag
  endif
endfor
h.photoz /= h.ngals

ii = where(h.out)
srt = reverse(sort(h(ii).ngals))
for i = 0L, N_ELEMENTS(ii) - 1 do begin
  h(ii(srt(i))).rank = i+1
endfor
gals(*).rank = -1
for i = 0L, ngout - 1 do begin
  ih = gals(i).haloID
  if (ih lt 0) then continue
  gals(i).rank = h(ih).rank
endfor

openw,1,'~/ki01/test_clusters.ascii.junk'
for i = 0L, N_ELEMENTS(ii) - 1 do begin
  ih = ii(srt(i))
  printf,1,h(ih).rank, h(ih).ra, h(ih).dec, h(ih).photoz
endfor
close,1

openw,1,'~/ki01/test_galaxies.ascii.junk'
for i = 0L, ngout - 1 do begin
;  if (gals(i).rank eq -1) then continue
  printf,1,gals(i).rank, gals(i).galaxyID, gals(i).photoz
endfor
close,1

end
