#!/usr/bin/env python
from __future__ import print_function, division
from rdel import luminosityfunction, config
from glob import glob
from copy import deepcopy as copy
import numpy as np
import shutil
import argparse
import pickle
import sys
import os


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('addgalsnames')
    parser.add_argument('cfgfile')
    parser.add_argument('--rdmodelfile')
    parser.add_argument('--evol')
    parser.add_argument('--Q')
    parser.add_argument('--z_redfraction1')
    parser.add_argument('--z_redfraction2')
    parser.add_argument('--redfraction1')
    parser.add_argument('--redfraction2')
    parser.add_argument('--scatter')
    parser.add_argument('--lbcgname')
    parser.add_argument('--hv')

    args = parser.parse_args()
    addgalsnames = args.addgalsnames
    cfgfile      = args.cfgfile

    cfg        = config.readCfg(cfgfile)
    print('parsing config')
    sim, lf, m = config.parseConfig(cfg)
    print('done')

    lbcgmodels = np.array(glob('{}/rdel/*lc*pkl'.format(cfg['Simulation']['outdir'])))
    scales     = np.array([float(l.split('_')[-2].split('.list')[0]) for l in lbcgmodels])
    zs         = 1/scales - 1
    
    zidx       = np.argsort(zs)
    lbcgmodels = lbcgmodels[zidx]
    zs         = zs[zidx]
    
    paramfiles = glob(addgalsnames)

    dt    = np.dtype([('param', 'S200'), ('value',float)])
    sdt    = np.dtype([('param', 'S200'), ('value','S200')])
    mbins = np.linspace(-25,0,2500)

    if args.rdmodelfile is not None:
        rdelparams = np.genfromtxt(args.rdmodelfile, dtype=dt)
    else:
        rdelparams = None

    for f in paramfiles:
        base = f.split('/')
        base = '/'.join(base[:-1])
        
        params = np.genfromtxt(f, dtype=dt)

#        if args.lbcgname is not None:
        sname = '{}/StringParameters'.format(base)
        sparams = np.genfromtxt(sname, dtype=sdt)
        sparams['value'][sparams['param']=='lbcgfile'] = '{}/lbcg_mvir_model.txt'.format(base)
        np.savetxt(sname, sparams, fmt='%s')
            
        zmin   = params['value'][params['param']=='ZREDMIN']
        zmax   = params['value'][params['param']=='ZREDMAX']

        if rdelparams is not None:
            psidx,   = np.where(params['param']=='Magmin')
            peidx,   = np.where(params['param']=='BCG_Mass_lim')

            psidx = psidx[0]
            peidx = peidx[0]
        
            params = np.hstack([params[:psidx+1], rdelparams, params[peidx:]])

            #params['value'][psidx] = (float(params['value'][psidx]) + 1.618 * (1/(1 + zmin) - 1/1.1))

        if args.Q is not None:
            params['value'][params['param']=='Q'] = args.Q
        if args.evol is not None:
            params['value'][params['param']=='evolution'] = args.evol
        if args.z_redfraction1 is not None:
            params['value'][params['param']=='Z_REDFRACTION1'] = args.z_redfraction1
        if args.z_redfraction2 is not None:
            params['value'][params['param']=='Z_REDFRACTION2'] = args.z_redfraction2
        if args.redfraction2 is not None:
            params['value'][params['param']=='REDFRACTION2'] = args.redfraction2
        if args.redfraction2 is not None:
            params['value'][params['param']=='REDFRACTION2'] = args.redfraction2
        if args.scatter is not None:
            params['value'][params['param']=='SCATTER'] = args.scatter

        np.savetxt(f, params, fmt='%s')
        
        zmean  = (zmax + zmin) / 2
        lfz    = lf.genLuminosityFunctionZ(mbins, zmean)
        lfz = lfz.view((np.float,2))

        if lf.unitmap['mag'] == 'mag':
            lfz[:,0] -= 5 * np.log10(sim[0].h)

        zidx = zs.searchsorted(zmean)

        print('Using model {}'.format(lbcgmodels[zidx]))

        with open(lbcgmodels[zidx], 'r') as fp:
            lbcgm = pickle.load(fp)

        print(lbcgm['lcmass_params'])
        pars = np.atleast_2d(lbcgm['lcmass_params'])
        np.savetxt('{}/lbcg_mvir_model.txt'.format(base), pars, fmt=['%6f','%6f','%6f','%6f','%6f'])
        np.savetxt('{}/LF.dat'.format(base), lfz, fmt=['%6f', '%6e'])

        if args.hv is not None:
            shutil.copyfile(args.hv, '{}/hv'.format(base))
