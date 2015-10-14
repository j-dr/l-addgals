pro get_millennium_hids, group_info_file, dbase, fbase, num, outfile

;base = '/nfs/slac/g/ki/ki01/mbusha/projects/ubercomp/BaseData/Simulation/'
;group_info_file = base+'analysis/group_info/dataoutput/group_info_006'
;outfile = base+'analysis/group_info/dataoutput/ahid_006'
;dbase = base+'output/'
;fbase = 'lcdm_fb.downsample.gad_'
;num = '006'

;readgroupinfo, group_info_file, ng, pos=gpos
rdfloat, group_info_file, m200, sigma, np, hx,hy,hz
m200 = 0
sigma = 0
np = 0
ng = N_ELEMENTS(hx)
gpos = fltarr(3,ng)
gpos(0,*) = hx
gpos(1,*) = hy
gpos(2,*) = hz
hx = 0
hy = 0
hz = 0
readgadget2, dbase, fbase, num, npart=npart, pos=pos,box=box
np = npart(1)
ahid = lonarr(np)
tgpos = fltarr(3,ng)
boxhalf = 0.5*box
print, "looping through particles..."
next = 0.0
for i = 0L, np - 1 do begin
  if (float(i)/np ge next) then begin
    percent = string(float(i)/np*100)
    print, percent+"% done"
    next += 0.01
  endif
  min_dist = 1e9
  for j = 0L, ng - 1 do begin
    tdx = gpos(0,j)-pos(0,i)
    if (tdx gt boxhalf) then tdx -= box
    if (tdx lt -boxhalf) then tdx += box
    if (abs(tdx) gt min_dist) then continue
    tdy = gpos(1,j)-pos(1,i)
    if (tdy gt boxhalf) then tdy -= box
    if (tdy lt -boxhalf) then tdy += box
    if (abs(tdy) gt min_dist) then continue
    tdz = gpos(2,j)-pos(2,i)
    if (tdz gt boxhalf) then tdz -= box
    if (tdz lt -boxhalf) then tdz += box
    if (abs(tdz) gt min_dist) then continue
    this_dist = sqrt(tdx*tdx + tdy*tdy + tdz*tdz)
    if (this_dist lt min_dist) then begin
      min_dist = this_dist
      ahid(i) = j
    endif
  endfor
endfor

openw,1,outfile
writeu,1,np
writeu,1,ahid
close,1
end
