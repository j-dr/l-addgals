pro run_zcarlos_dr8
;pro run_zcarlos_dr8, pixel

args=command_line_args(count=ct)
pixel = args[0]
;gfile = '../../truth/Aardvark_v0.5_truth.'+pixel+'.fit'
;dr8_training_set = '../../photoz_DR8/PO_Aardvark_DR8_training_set_sdss_mag.fit'
;dr8_obsfile = '../../DR8/Aardvark_v0.5_sdss_mag.'+pixel+'.fit'
;dr8_ofile = '../../photoz_DR8/Aardvark_v0.5_DR8_zCarlos.'+pixel+'.fit'
gfile = args[1]
dr8_training_set = args[2]
dr8_obsfile = args[3]
dr8_ofile = args[4]
;des_training_set = '../../truth/BCC_Optimistic_training_set.fit'
;des_obsfile = '../../truth/Aardvark_v0.5_truth.'+pixel+'.fit'
;des_ofile = '../../photoz/Aardvark_v0.5_zCarlos.'+pixel+'.fit'
;des_training_set = args[5]
;des_obsfile = args[6]
;des_ofile = args[7]

;;;first we do the DR8 catalog
print, gfile
g = mrdfits(gfile,1)
make_carlos_photoz, dr8_training_set, dr8_obsfile, maglim = 21.8, $
        limband = 2, zmin = 0.0, zmax = 1.1, $
        grid = '35', res = '35', ofile=dr8_ofile, id=g.id, $
        ophotnum5 = 'photnum_DR8_nne5.'+pixel+'.num', $
        ophotnum100 = 'photnum_DR8_nne100.'+pixel+'.num', $
        onnweight5 = 'nnweight_DR8_nne5.'+pixel+'.prw', $
        onnweight100 = 'nnweight_DR8_nne100.'+pixel+'.prw'


;;;second we do the DES optimistic catalog
;make_carlos_photoz, des_training_set, des_obsfile, maglim = 24.0, $
;        limband = 2, zmin = 0.0, zmax = 3.0, $
;        grid = '100', res = '100', ofile=des_ofile, $
;        ophotnum5 = 'photnum_DES_nne5.'+pixel+'.num', $
;        ophotnum100 = 'photnum_DES_nne100.'+pixel+'.num', $
;        onnweight5 = 'nnweight_DES_nne5.'+pixel+'.prw', $
;        onnweight100 = 'nnweight_DES_nne100.'+pixel+'.prw'



end

