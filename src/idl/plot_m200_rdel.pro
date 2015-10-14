@~/mbusha.1/lib/RisaLibs/idl/mylib/plothist

base = '/nfs/slac/g/ki/ki02/mbusha/projects/addgals/Millennium/downsample100/analysis/'

rnn_file12 = base+'rnn/dataoutput/rnn_downsample100_063_nn21'
rnn_file13 = base+'rnn/dataoutput/rnn_downsample100_063_nn209'
ahid_file = base+'group_info/dataoutput/ahid_downsample100_063.rnn'
halo_file = base+'group_info/dataoutput/halos_063_all'

rhocrit = 3*100*100/(8*!PI*4.301e-9)

;readrnn, rnn_file13, np=np13, rnn=rnn13
;readrnn, rnn_file12, np=np12, rnn=rnn12
;rdfloat, ahid_file, ahid, rhalo
;ahid = long(ahid)
;rdfloat, halo_file, m200
r200 = (m200/(200.*rhocrit*4./3*!PI))^(1./3)

np = np12
p_m200 = fltarr(np)
p_r200 = fltarr(np)
p_m200(*) = m200(ahid)
p_r200(*) = r200(ahid)
out = where(rhalo gt p_r200, comp = in)
;p_m200(out) = 1e8
r12 = (1.8e12/(4./3*!PI*200*rhocrit))^(1./3)
r13 = (1.8e13/(4./3*!PI*200*rhocrit))^(1./3)

ran = floor(randomu(2719,100000)*N_ELEMENTS(in))
window, 0
plot, rnn13(in(ran)), p_m200(in(ran)), psym=3, /ylog, /ysty, $
      xtitle = 'R!Ddelta', ytitle = 'M!D200'
oplot, rnn12(in(ran)), p_m200(in(ran)), psym=3, color=!red
oplot, [0,10], 1.8e13*[1,1]
oplot, [0,10], 1.8e12*[1,1], color=!red
oplot, r13*[1,1], [1,1e18]
oplot, r12*[1,1], [1,1e18], color=!red

window,1
plothist, rnn13(out), bin=0.1, yrange = [0,6e6], xtitle = 'R!Ddelta'
plothist, rnn13(in), bin=0.1, /over, linestyle=1
plothist, rnn12(out), bin=0.1, color=!red, /over
plothist, rnn12(in), bin=0.1, /over, linestyle=1, color=!red
oplot, r13*[1,1], [0,1e12], linestyle = 2
oplot, r12*[1,1], [0,1e12], linestyle = 2, color=!red

end
