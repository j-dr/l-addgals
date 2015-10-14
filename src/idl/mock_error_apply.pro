pro mock_error_apply,mode,tmag,oflux,ofluxerr,omag,omagerr,seed=seed, point_source=point_source

if n_params() eq 0 then begin
    print,'syntax- mock_error_apply,mode,tmag,oflux,ofluxerr,omag,omagerr,seed=seed'
    print,' mode is DR8, STRIPE82, CFHTLS, DEEP2, FLAMEX, IRAC, NDWFS, RCS'
    print,' takes an array of tmags, returns arrays of omags, omagerrs'
    return
endif

if n_elements(seed) eq 0 then seed=systime(/seconds)


case (strupcase(mode)) of
    'DR8':begin
        ;; u,g,r,i,z
        maglims = [20.425,21.749,21.239,20.769,19.344]
        exptimes = [21.00,159.00,126.00,99.00,15.00]
        lnscat = [0.284,0.241,0.229,0.251,0.264]
    end
    'STRIPE82':begin
        ;; u,g,r,i,z
        maglims = [22.070,23.432,23.095,22.649,21.160]
        exptimes = [99.00,1172.00,1028.00,665.00,138.00]
        lnscat = [0.280,0.229,0.202,0.204,0.246]
    end
    'CFHTLS':begin
        ;; u,g,r,i,z
        ;; use W1, MAG_AUTO
        maglims = [24.298,24.667,24.010,23.702,22.568]
        exptimes = [2866.00,7003.00,4108.00,3777.00,885.00]
        lnscat = [0.259,0.244,0.282,0.258,0.273]
    end
    'DEEP2':begin
        ;; B,R,I
        maglims = [24.730,24.623,24.092]
        exptimes = [7677.00,8979.00,4402.00]
        lnscat = [0.300,0.293,0.300]
    end
    'FLAMEX':begin
        ;; J, K_s
        maglims = [21.234,20.929]
        exptimes = [259.00,135.00]
        lnscat = [0.300,0.289]
    end
    'IRAC':begin
        ;; 3.4 um, 4.6 um
        maglims = [19.352,18.574]
        exptimes = [8.54,3.46]
        lnscat = [0.214,0.283]
    end
    'NDWFS':begin
        ;; B_w, R, I
        maglims = [25.142,23.761,23.650]
        exptimes = [6391.00,1985.00,1617.00]
        lnscat = [0.140,0.294,0.272]
    end
    'RCS':begin
        ;; g, r, i, z
        maglims = [23.939,23.826,23.067,21.889]
        exptimes = [2850.00,2568.00,1277.00,431.00]
        lnscat = [0.164,0.222,0.250,0.271]
    end
    'VHS':begin
        ;; J, H, K_s
        maglims = [20.141,19.732,19.478]
        exptimes = [36.00,31.00,23.00]
        lnscat = [0.097,0.059,0.069]
    end
    'VIKING':begin
        ;; z, y, J, H, K_s
        maglims = [21.643,20.915,20.768,20.243,20.227]
        exptimes = [622.00,246.00,383.00,238.00,213.00]
        lnscat = [0.034,0.048,0.052,0.040,0.066]
    end
    'DC6B':begin
        ;; g,r,i,z  (no Y)
        maglims = [24.486,23.473,22.761,22.402]
        exptimes = [2379.00,1169.00,806.00,639.00]
        lnscat = [0.300,0.300,0.300,0.300]
    end
    'DES': begin
        maglims = [24.956,24.453,23.751,23.249,21.459]
        exptimes = [14467.00,12471.00,6296.00,5362.00,728.00]
        lnscat = [0.2,0.2,0.2,0.2,0.2]
    end
    'BCS_LO':begin
        maglims = [23.082,22.618,22.500,21.065]
        exptimes = [809.00,844.00,641.00,108.00]
        lnscat = [0.277,0.284,0.284,0.300]
    end
    'BCS':begin
        maglims = [23.360,23.117,22.539,21.335]
        exptimes = [838.00,1252.00,772.00,98.00]
        lnscat = [0.276,0.272,0.278,0.279]
    end
    'DES_SV': begin
	maglims = [23.621,23.232,23.008,22.374,20.663]
	exptimes = [4389.00,1329.00,1405.00,517.00,460.00]
	lnscat = [0.276,0.257,0.247,0.241,0.300] 
    end
    'DES_SV_OPTIMISTIC': begin
        maglims = [23.621+0.5,23.232+0.5,23.008,22.374,20.663]
        exptimes = [4389.00,1329.00,1405.00,517.00,460.00]
        lnscat = [0.276,0.257,0.247,0.241,0.300]
    end
    'WISE': begin
        ;; 3.4 um, 4.6 um
        maglims = [19.352,18.574]
        exptimes = [8.54,3.46]
        lnscat = [0.214,0.283]
    end
    'DECALS': begin
       maglims = [23.3,23.3,22.2,20.6,19.9]
       exptimes = [1000,3000,2000,1500,1500]
       lnscat = [0.2,0.2,0.2,0.2,0.2]
    end
    else: begin
        print,'Illegal mode: ',mode
        return
    end
endcase

nmag=n_elements(maglims)

if (nmag ne n_elements(tmag[*,0])) then begin
    print,'Error: ',mode,' requires ',nmag,' magnitudes, but input tmag array has ',n_elements(tmag[*,0]),' magnitudes.'
    return
endif

ngal=n_elements(tmag[0,*])

zp=22.5

;; calculate fsky1 -- sky in 1 second
print, 'flux1_lim'
flux1_lim = 10.^((maglims-zp)/(-2.5)) > 120/exptimes ;; don't go negative...
print, 'fsky1'
fsky1 = (flux1_lim^2.*exptimes)/100. - flux1_lim

print, 'ostuff'
oflux=fltarr(nmag,ngal)
ofluxerr=oflux
omag=oflux
omagerr=oflux

print, 'offset'
offset = 0.0
if KEYWORD_SET(point_source) then begin
    offset = 0.6
    lnscat /= 4.0
endif

print, 'nmag'
for i=0l,nmag-1 do begin
    tflux = exptimes[i] * 10.^((reform(tmag[i,*])-offset-zp)/(-2.5))

    noise = exp(alog(sqrt(fsky1[i]*exptimes[i] + tflux)) + lnscat[i]*randomn(seed,ngal))

    flux = tflux + noise*randomn(seed,ngal)

    oflux[i,*] = flux / exptimes[i]
    ofluxerr[i,*] = noise/exptimes[i]

    omag[i,*] = 22.5-2.5*alog10(oflux[i,*])
    omagerr[i,*] = (2.5/alog(10.))*(ofluxerr[i,*]/oflux[i,*])

    bad=where(finite(omag[i,*]) eq 0,nbad)
    if (nbad gt 0) then begin
        omag[i,bad] = 99.0
        omagerr[i,bad] = 99.0
    endif

endfor


return
end


















