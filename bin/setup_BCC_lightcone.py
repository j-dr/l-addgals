#!/usr/bin/env python
from __future__ import print_function, division
from rdel import luminosityfunction, config
from glob import glob
import numpy as np
import pickle
import sys
import os


if __name__ == '__main__':

    addgalsnames = sys.argv[1]
    cfgfile     = sys.argv[2]

    cfg        = config.readCfg(cfgfile)
    sim, lf, m = config.parseConfig(cfg)

    lbcgmodels = np.array(glob('{}/rdel/*lc*pkl'.format(cfg['Simulation']['outdir'])))
    scales     = np.array([float(l.split('_')[-2].split('.list')[0]) for l in lbcgmodels])
    zs         = 1/scales - 1
    
    zidx       = np.argsort(zs)
    lbcgmodels = lbcgmodels[zidx]
    zs         = zs[zidx]
    
    paramfiles = glob(addgalsnames)

    dt    = np.dtype([('param', 'S20'), ('value',float)])
    mbins = np.linspace(-25,0,2500)

    for f in paramfiles:
        base = f.split('/')
        base = '/'.join(base[:-1])
        
        params = np.genfromtxt(f, dtype=dt)
        zmin   = params['value'][params['param']=='ZREDMIN']
        zmax   = params['value'][params['param']=='ZREDMAX']

        zmean  = (zmax + zmin) / 2
        lfz    = lf.genLuminosityFunctionZ(mbins, zmean)
        lfz = lfz.view((np.float,2))

        zidx = zs.searchsorted(zmean)
        with open(lbcgmodels[zidx], 'r') as fp:
            lbcgm = pickle.load(fp)

        print(lbcgm['lcmass_params'])
        print(lbcgm['lcmass_params'].shape)
        pars = np.atleast_2d(lbcgm['lcmass_params'])
        np.savetxt('{}/lbcg_mvir_model.txt'.format(base), pars, fmt=['%6f','%6f','%6f','%6f','%6f'])
        np.savetxt('{}/LF.dat'.format(base), lfz, fmt=['%6f', '%6e'])

        
