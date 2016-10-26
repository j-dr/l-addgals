from __future__ import print_function, division
from mpi4py import MPI
from glob import glob
from copy import copy
from itertools import izip
import numpy as np
import healpy as hp
import healpix_util as hu
import fitsio
import os
import sys
import yaml

MASKED_VAL = -9999

def read_partial_map(filename, ext=1, masked_val=MASKED_VAL,
                        pix_col='PIXEL', val_col='SIGNAL'):
    f=fitsio.FITS(filename)[ext]
    nside=int(f.read_header()['NSIDE'])
    hpix=hu.HealPix("ring",nside)
    m=masked_val*np.ones(hpix.npix)
    m_data=f.read()
    m[m_data[pix_col]]=m_data[val_col]
    return hu.Map("ring", m)

def sys_map_cuts(gal_data, sys_map_data=None, ra_col='ra',
                  dec_col='dec'):
    """Get systematic values, and cut data without sys data"""

    sys_map_vals={}
    mask=np.ones(len(gal_data),dtype='bool')

    for i, (name, m) in enumerate(sys_map_data.iteritems()):

         sys_map_vals[name]=m.get_mapval(gal_data[ra_col],
                                            gal_data[dec_col])

         if i==0:
             use = sys_map_vals[name]!=MASKED_VAL
         else:
             use &= sys_map_vals[name]!=MASKED_VAL

    return use, sys_map_vals

def gold_cuts(gal_data, ra_col='RA', dec_col='DEC',
                gold_fp_map=None, gold_br_map=None):
    if gold_fp_map is None:
        gold_fp_map=hu.readMap(gold_footprint_fn)
    if gold_br_map is None:
        gold_br_map=hu.readMap(gold_badreg_fn)
    use=((gold_fp_map.get_mapval(gal_data[ra_col],gal_data[dec_col])>=1) *
             (gold_br_map.get_mapval(gal_data[ra_col],gal_data[dec_col])==0)) #LSS uses < 3?
    return use.astype(bool)


def WL_cuts(obs, truth, pz, sys_map_vals,
                maglim_cut_factor, rgrp_cut, z_col):

    psf_size = 0.26 * 0.5 * sys_map_vals['psf_fwhm_r']
    mag_mask = obs['MAG_R'] < -2.5*np.log10(maglim_cut_factor) + sys_map_vals['maglim_r']
    size_mask = np.sqrt(obs['SIZE']**2 + psf_size**2) > rgrp_cut * psf_size

    if z_col=='Z':
        z = truth[z_col]
    else:
        z = pz[z_col]

    other_mask = (np.isfinite(sys_map_vals['maglim_r']) * np.isfinite(sys_map_vals['psf_fwhm_r'])
                  * (z > 0.1))

    good = mag_mask * size_mask * other_mask

    return good

def LSS_cuts(obs, truth, pz, sys_map_vals, zcol):

    if 'MEAN_Z' == zcol:
        z = pz[zcol]
    else:
        z = truth[zcol]

    mask = (obs['MAG_I'] > 17.5) & \
            (obs['MAG_I'] < 22)  & \
            (obs['MAG_I'] < (19.0 + 3.0 * z)) & \
            ((obs['MAG_I'] - obs['MAG_Z'] + 2.0 * obs['MAG_R'] - obs['MAG_I']) > 1.7) & \
            (-1 < (obs['MAG_G'] - obs['MAG_R'])) & \
            ((obs['MAG_G'] - obs['MAG_R']) < 3)   & \
            (-1 < (obs['MAG_R'] - obs['MAG_I'])) & \
            ((obs['MAG_R'] - obs['MAG_I']) < 2.5)   & \
            (-1 < (obs['MAG_I'] - obs['MAG_Z'])) & \
            ((obs['MAG_I'] - obs['MAG_Z']) < 2.)   & \
            ((obs['RA'] < 15.) | (obs['RA']>290) | (obs['DEC']<-35))

    return mask

def make_single_selection(obs, truth, pz, mask,
                            sys_map_data, cut_func, zcol):

    #mask based on systematics maps
    print('Survey mask: {}'.format(mask))
    print('{}'.format(mask.any()))
    smask, sys_map_vals = sys_map_cuts(obs,
                            sys_map_data=sys_map_data,
                            ra_col='RA',dec_col='DEC')

    print('Systematic mask: {}'.format(smask))
    print('{}'.format(smask.any()))    
    mask &= smask

    #apply galaxy property cuts
    mask &= cut_func(obs, truth, pz, sys_map_vals, zcol)
    print('Cut mask: {}'.format(mask))
    print('{}'.format(mask.any()))
    
    return mask

def pair_files(ofiles, tfiles, pzfiles):

    opix = np.array([int(i.split('.')[-2]) for i in ofiles])
    tpix = np.array([int(i.split('.')[-2]) for i in tfiles])
    ppix = np.array([int(i.split('.')[-2]) for i in pzfiles])

    ssidx = np.in1d(tpix, opix)

    tfiles = tfiles[ssidx]
    tpix   = tpix[ssidx]

    assert(len(tpix)==len(opix))
    assert(len(tpix)==len(ppix))

    oidx = opix.argsort()
    tidx = tpix.argsort()
    pidx = ppix.argsort()

    ofiles = ofiles[oidx]
    tfiles = tfiles[tidx]
    pzfiles = pzfiles[pidx]

    return ofiles, tfiles, pzfiles


if __name__=="__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    cfgfile = sys.argv[1]
    with open(cfgfile, 'r') as fp:
        cfg = yaml.load(fp)
    print(cfg)
    #Read in gold masks
    gold_fp=hu.readMap(cfg['gold']['gold_footprint_fn'])
    gold_br=hu.readMap(cfg['gold']['gold_badreg_fn'])

    ofiles = np.array(glob(cfg['sim']['obspath']))
    tfiles = np.array(glob(cfg['sim']['truthpath']))
    pzfiles = np.array(glob(cfg['sim']['pzpath']))

    ofiles, tfiles, pzfiles = pair_files(ofiles, tfiles, pzfiles)

    sys_map_data = {}

    for of, tf, pz in izip(ofiles[rank::size], tfiles[rank::size], pzfiles[rank::size]):

        obsf   = fitsio.FITS(of, 'rw')
        truthf = fitsio.FITS(tf, 'r')
        pf     = fitsio.FITS(pz, 'r')

        obs    = obsf[-1].read()
        truth  = truthf[-1].read()
        pz     = pf[-1].read()

        mask = gold_cuts(obs, gold_fp_map=gold_fp,
                            gold_br_map=gold_br)

        for sample in cfg['samples']:
            smask = copy(mask)
            scfg = cfg['samples'][sample]
            if 'sys_maps' in scfg.keys():
                for name, mfile in scfg['sys_maps'].iteritems():
                    if name not in sys_map_data.keys():
                        m=read_partial_map(mfile,masked_val=np.nan)
                        sys_map_data[name] = m

            if sample == 'LSS':
                cut_fcn = LSS_cuts
            elif sample == 'WL':
                cut_fcn = lambda o, t, p, sys_map_vals, zcol : WL_cuts(o, t, p, 
                                                                   sys_map_vals,
                    scfg['maglim_cut_factor'],
                    scfg['rgrp_cut'], zcol)

            print(sys_map_data)

            smask = make_single_selection(obs, truth, pz, mask, sys_map_data, cut_fcn, scfg['z_col'])
            print('Any in {} sample?: {}'.format(sample, smask.any()))
            #add flag to file here!
            sflag = '{0}_FLAG'.format(sample)
            if sflag not in obs.dtype.names:
                obsf[-1].insert_column(sflag, smask.astype(int))
            else:
                obsf[-1].write_column(sflag, smask.astype(int))

        obsf.close()
        truthf.close()
        pf.close()
