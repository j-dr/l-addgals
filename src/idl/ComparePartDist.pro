;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PROGRAM: ComparePartDist
; PURPOSE: Looks at distances from particles to their linked halos
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

HaloFile = '/nfs/slac/g/ki/ki01/mbusha/data/hubble/lcdm.MS.msort12'
rdfloat, HaloFile, m200, zred, sigma, ip, x, y, z, vx, vy, vz, siglos, rdelta, id, skipline=2
m200 *= 1e15
x *= 1000.
y *= 1000.
z *= 1000.

pbase = '/nfs/slac/g/ki/ki01/mbusha/data/hubble/cubes/position/p.'
ahidbase = '/nfs/slac/g/ki/ki01/mbusha/data/hubble/rnn/ahid.'
for i = 8, 11 do begin
  num = '00'+strcompress(string(i), /remove_all)
  num = strmid(num, strlen(num)-2,2)
  pfile = pbase+num+'.'+num+'.'+num
  ahidfile = ahidbase+num+'.'+num+'.'+num
  rdfloat, ahidfile, tahid
  np = N_ELEMENTS(tahid)
  tind = intarr(3)
  tpos = fltarr(3,np)
  openr,1,pfile
  for i = 0L, np - 1 do begin
      readu,1,tind
      tpos(*,i) = tind(*)
