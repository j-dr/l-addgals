pro run_zcarlos_bossanova

args=command_line_args(count=ct)
pixel = args[0]
gfile = args[1]
dr8_training_set = args[2]
dr8_obsfile = args[3]
dr8_ofile = args[4]
label = args[5]

print, gfile
g = mrdfits(gfile,1)
make_carlos_photoz, dr8_training_set, dr8_obsfile, maglim = 21.8, $
        limband = 2, zmin = 0.0, zmax = 1.1, $
        grid = '35', res = '35', ofile=dr8_ofile, id=g.id, $
        ophotnum5 = 'photnum_'+label+'_nne5.'+pixel+'.num', $
        ophotnum100 = 'photnum_'+label+'_nne100.'+pixel+'.num', $
        onnweight5 = 'nnweight_'+label+'_nne5.'+pixel+'.prw', $
        onnweight100 = 'nnweight_'+label+'_nne100.'+pixel+'.prw'

end

