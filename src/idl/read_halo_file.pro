function read_halo_file, hfile

command = 'wc '+hfile
spawn, command, result
nint = strsplit(result)
nhalos = long(strmid(result,0,nint(1)))

s1=create_struct('haloid',0L,'z',0.0,'m200',0.0, 'r200',0.0, 'ngals',0, $
                 'ra',0.0,'dec',0.0,'halopx',0.0,'halopy',0.0,'halopz',0.0,$
                 'halovx',0.0,'halovy',0.0,'halovz',0.0,'siglos',0.0)
halos=replicate(s1,nhalos)

openr,1,hfile
for i = 0L, nhalos(0) - 1 do begin
  readf,1,t1,t2,t3,t4,t5,t6,t7
  halos(i).haloid = t1
  halos(i).m200 = t2
  halos(i).r200 = t3
  halos(i).ra = t4
  halos(i).dec = t5
  halos(i).z = t6
  halos(i).ngals = t7
endfor
close,1

HV_file = 'HV_halos.dat'
command = 'wc '+HV_file
spawn, command, result
nint = strsplit(result)
nHVhalos = long(strmid(result,0,nint(1))) - 2
s2 = create_struct('m15',0.0,'zred',0.0,'sigma',0,'ip',0,'xGpc',0.0,$
                   'yGpc',0.0,'zGpc',0.0,'vx',0,'vy',0,'vz',0,'siglos',0,$
                   'rdelta',0.0,'id',0L)
HVhalos = replicate(s2,nHVhalos)
junk = 'asdf'
openr,1,HV_file
readf,1,junk
readf,1,junk
for i = 0L, nHVhalos(0) - 1 do begin
  readf,1,s2
  HVhalos(i) = s2
endfor
for i = 0L, nhalos(0) - 1 do begin
  th = halos(i).haloid-1
  halos(i).halopx = HVhalos(th).xGpc*1000.
  halos(i).halopy = HVhalos(th).yGpc*1000.
  halos(i).halopz = HVhalos(th).zGpc*1000.
  halos(i).halovx = HVhalos(th).vx
  halos(i).halovy = HVhalos(th).vy
  halos(i).halovz = HVhalos(th).vz
  halos(i).siglos = HVhalos(th).siglos
endfor
return, halos

END
