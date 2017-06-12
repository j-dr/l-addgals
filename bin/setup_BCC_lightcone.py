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
    parser.add_argument('--QBASELINE')
    parser.add_argument('--EVOLVECEN')
    parser.add_argument('--rfmodelfile')
    parser.add_argument('--scatter')
    parser.add_argument('--lbcgname')
    parser.add_argument('--hv')
    parser.add_argument('--Magmin_dens')

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
    print(lbcgmodels)
    zs         = zs[zidx]
    
    paramfiles = glob(addgalsnames)

    dt    = np.dtype([('param', 'S200'), ('value',float)])
    sdt    = np.dtype([('param', 'S200'), ('value','S200')])
    mbins = np.linspace(-25,0,2500)

    if args.rdmodelfile is not None:
        rdelparams = np.genfromtxt(args.rdmodelfile, dtype=dt)
    else:
        rdelparams = None

    if args.rfmodelfile is not None:
        rfparams = np.genfromtxt(args.rfmodelfile, dtype=dt)
    else:
        rfparams = None

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

        if rfparams is not None:
            print(rfparams['param'])
            if rfparams['param'][0] in params['param']:
                print(params['param']==rfparams['value'][0])
                print(rfparams['value'][0])
                psidx, = np.where(params['param']==rfparams['param'][0])
                peidx, = np.where(params['param']==rfparams['param'][-1])
                psidx  = psidx[0]
                peidx  = peidx[-1]
                params = np.hstack([params[:psidx], rfparams, params[peidx+1:]])
            elif 'REDFRACTION1' in params['param']:
                psidx, = np.where(params['param']=='REDFRACTION1')
                peidx, = np.where(params['param']=='Z_REDFRACTION2')
                psidx  = psidx[0]
                peidx  = peidx[0]
                params = np.hstack([params[:psidx], rfparams, params[peidx+1:]])
            else:
                params = np.hstack([params, rfparams])
            
        if args.QBASELINE is not None:
            if 'QBASELINE' in params['param']:
                params['value'][params['param']=='QBASELINE'] = args.QBASELINE
            else:
                qbarr = np.zeros(1, dtype=dt)
                qbarr['param'][0] = 'QBASELINE'
                qbarr['value'][0] = args.QBASELINE
                params = np.hstack([params, qbarr])

        if args.Magmin_dens is not None:
            if 'Magmin_dens' in params['param']:
                params['value'][params['param']=='Magmin_dens'] = args.Magmin_dens
            else:
                qbarr = np.zeros(1, dtype=dt)
                qbarr['param'][0] = 'Magmin_dens'
                qbarr['value'][0] = args.Magmin_dens
                params = np.hstack([params, qbarr])
                
        if args.EVOLVECEN is not None:
            if 'EVOLVECEN' in params['param']:
                params['value'][params['param']=='EVOLVECEN'] = args.EVOLVECEN
            else:
                qbarr = np.zeros(1, dtype=dt)
                qbarr['param'][0] = 'EVOLVECEN'
                qbarr['value'][0] = args.EVOLVECEN
                params = np.hstack([params, qbarr])
                
        if args.Q is not None:
            params['value'][params['param']=='Q'] = args.Q
        if args.evol is not None:
            params['value'][params['param']=='evolution'] = args.evol
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
