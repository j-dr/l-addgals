from __future__ import print_function, division
from mpi4py import MPI
from glob import glob
import numpy.lib.recfunctions as rf
import healpy as hp
import numpy as np
import fitsio
import yaml
import sys

def compute_lensing(g, shear, halos=False):


    lensfields = ['GAMMA1', 'GAMMA2', 'KAPPA', 'W', 'MU', 'TRA', 'TDEC', 'RA', 'DEC']

    if not halos:
        lensfields.append("LMAG")

    nmimg  = 0
    adtype = []
    fields = []
    data   = []
    
    for f in lensfields:
        if f not in g.dtype.names:
            fields.append(f)
            print(f)
            if f == 'RA':
                adtype.append(np.dtype([(f,np.float)]))
                theta, phi = hp.vec2ang(g[['PX','PY','PZ']].view((g.dtype['PX'],3)))
                dec, ra =  -np.degrees(theta-np.pi/2.), np.degrees(np.pi*2.-phi)
                data.append(ra)
            elif f == 'DEC':
                adtype.append(np.dtype([(f,np.float)]))
                theta, phi = hp.vec2ang(g[['PX','PY','PZ']].view((g.dtype['PX'],3)))
                dec, ra =  -np.degrees(theta-np.pi/2.), np.degrees(np.pi*2.-phi)
                data.append(dec)
            elif f == 'LMAG':
                adtype.append(np.dtype([(f,(np.float,5))]))
                data.append(np.zeros(len(g)))
            else:
                adtype.append(np.dtype([(f,np.float)]))
                data.append(np.zeros(len(g)))

    if len(fields)>0:
        g = rf.append_fields(g, fields, data=data,
                             dtypes=adtype, usemask=False)

    print("Number of rows in shear catalog: {0}".format(len(shear)))
    print("Number of rows in galaxy catalog: {0}".format(len(g)))

    for i in range(len(shear)):
        tind = shear['index'][i]
        
        if (g[tind]['GAMMA1'] != -1):
            nmimg += 1
            g = np.hstack([g, g[tind]])
            tind = len(g)-1

        g['TRA'][tind] = g[tind]['RA']
        g['TDEC'][tind] = g[tind]['DEC']
        
        g['RA'][tind] = shear['ra'][i]
        g['DEC'][tind] = shear['dec'][i]
  
        #extract g1,g2,kappa,w from A using formulas from Vale & White (2003)
        # A = |1-kappa-gamma1      -gamma2-w|
        #     |-gamma2+w      1-kappa+gamma1|


        g['GAMMA1'][tind] = 0.5*(shear[i]['A11'] - shear[i]['A00'])
        g['GAMMA2'][tind] = -0.5*(shear[i]['A01'] + shear[i]['A10'])
        g['KAPPA'][tind] = 1.0 - 0.5*(shear[i]['A00'] + shear[i]['A11'])
        g['W'][tind] = 0.5*(shear[i]['A10'] - shear[i]['A01'])
    
        #compute mu = 1/detA
        g['MU'][tind] = 1./(shear[i]['A11']*shear[i]['A00'] - shear[i]['A01']*shear[i]['A10'])

        if not halos:
            #lens the size and magnitudes
            g['SIZE'][tind] = g[tind]['TSIZE']*np.sqrt(g[tind]['MU'])
            for im in range(g['AMAG'].shape[1]):
                g['LMAG'][tind,im] = g['TMAG'][tind,im] - 2.5*np.log10(g[tind]['MU'])
  
            #get intrinsic shape
            epss = complex(g['TE'][tind,0], g['TE'][tind,1])
  
            #;;get reduced complex shear g = (gamma1 + i*gamma2)/(1-kappa)
            gquot = complex(g['GAMMA1'][tind], g['GAMMA2'][tind]) / complex(1.-g['KAPPA'][tind], 0.)
  
            #;;lens the shapes - see Bartelmann & Schneider (2001), Section 4.2
            if (abs(gquot) < 1):
                eps = (epss+gquot)/(complex(1.0,0.0)+(epss*gquot.conjugate()))
            else:
                eps = (complex(1.0, 0.0)+(gquot*epss.conjugate()))/(epss.conjugate()+gquot.conjugate())

            g['EPSILON'][tind,0] = eps.real
            g['EPSILON'][tind,1] = eps.imag

    print("Number of multiply imaged galaxies: {0}".format(nmimg))

    return g


def add_lensing(gfiles, sfiles):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    assert(len(gfiles)==len(sfiles))

    gpix = np.array([int(gf.split('.')[-2]) for gf in gfiles])
    spix = np.array([int(sf.split('.')[-3].split('_')[0]) for sf in sfiles])

    gidx = gpix.argsort()
    sidx = spix.argsort()
    gpix = gpix[gidx]
    spix = spix[sidx]
    
    assert(all(gpix==spix))
    
    gfiles = gfiles[gidx]
    sfiles = sfiles[sidx]

    for i, f in enumerate(gfiles):
        if ('halo' in f) and ('halo' not in sfiles[i]):
            if ((i+1)<len(gpix)) and (gpix[i]==spix[i+1]):
                temp = sfiles[i+1]
                sfiles[i+1] = sfiles[i]
                sfiles[i] = temp
            elif gpix[i]==spix[i-1]:
                temp = sfiles[i-1]
                sfiles[i-1] = sfiles[i]
                sfiles[i] = temp

            assert(('halo' in f) and ('halo' in sfiles[i]))
                

    for gf,sf in zip(gfiles[rank::size],sfiles[rank::size]):

        print("[{1}] Lensing {0}".format(gf, rank))
        print("[{1}] Using  {0}".format(sf, rank))
        
        gfs = gf.split('/')
        gbase = "/".join(gfs[:-1])

        
        g     = fitsio.read(gf)
        shear = fitsio.read(sf)

        if 'halo' in gf:
            g     = compute_lensing(g, shear, halos=True)
        else:
            g     = compute_lensing(g, shear)

        gfss = gfs[-1].split('.')
        oname = "{0}/{1}_lensed.{2}.fits".format(gbase, gfss[0], gfss[1])

        print("[{1}] Writing lensed catalog to {0}".format(oname, rank))
        
        fits = fitsio.FITS(oname, 'rw')
        fits.write(g)


if __name__=='__main__':

    cfgfile = sys.argv[1]

    with open(cfgfile, 'r') as fp:
        cfg = yaml.load(fp)

    snames = np.loadtxt(cfg['LensGalsList'], dtype=str)
    gnames = np.loadtxt(cfg['TruthGalsList'], dtype=str)
  
    add_lensing(gnames, snames)
