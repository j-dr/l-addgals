from __future__ import print_function, division
from mpi4py import MPI
from glob import glob
import numpy as np
import numpy.lib.recfunctions as rf
import fitsio
import sys

def compute_lensing(g, shear):
  
    if 'W' not in g.dtype:
        adtype = [np.dtype([('W',np.float)])]
        data = [np.zeros(len(g))]
        g = rf.append_fields(g,['W'], data=data,
                             dtypes=adtype, usemask=False)

    for i in range(len(shear)):
        tind = shear['index'][i]

        if (g[tind]['GAMMA1'] != 0.):
            g = np.hstack([g, g[tind]])
            tind = len(g)-1

        g[tind].ra = shear[i].ra
        g[tind].dec = shear[i].dec
  
        #extract g1,g2,kappa,w from A using formulas from Vale & White (2003)
        # A = |1-kappa-gamma1      -gamma2-w|
        #     |-gamma2+w      1-kappa+gamma1|

        g[tind]['GAMMA1'] = 0.5*(shear[i]['A11'] - shear[i]['A00'])
        g[tind]['GAMMA2'] = -0.5*(shear[i]['A01'] + shear[i]['A10'])
        g[tind]['KAPPA'] = 1.0 - 0.5*(shear[i]['A00'] + shear[i]['A11'])
        g[tind]['W'] = 0.5*(shear[i]['A10'] - shear[i]['A01'])
    
        #compute mu = 1/detA
        g[tind]['MU'] = 1./(shear[i]['A11']*shear[i]['A00'] - shear[i]['A01']*shear[i]['A10'])
  
        #lens the size and magnitudes
        g[tind]['SIZE'] = g[tind]['TSIZE']*np.sqrt(g[tind]['MU'])
        for im in range(g['AMAG'].shape[1]):
            g[tind]['LMAG'][im] = g[tind]['TMAG'][im] - 2.5*np.log10(g[tind]['MU'])
  
        #get intrinsic shape
        epss = complex(g[tind]['TE'][0], g[tind]['TE'][1])
  
        #;;get reduced complex shear g = (gamma1 + i*gamma2)/(1-kappa)
        gquot = complex(g[tind]['GAMMA1'], g[tind]['GAMMA2']) / complex(1.-g[tind]['KAPPA'], 0.)
  
        #;;lens the shapes - see Bartelmann & Schneider (2001), Section 4.2
        if (abs(gquot) < 1):
            eps = (epss+gquot)/(complex(1.0,0.0)+(epss*gquot.conjugate()))
        else:
            eps = (complex(1.0, 0.0)+(gquot*epss.conjugate()))/(epss.conjugate()+gquot.conjugate())

        g[tind]['EPSILON'][0] = eps.real
        g[tind]['EPSILON'][1] = eps.imag

    return g


def add_lensing(gfiles, sbase):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    for gf in gfiles[rank::size]:
        print("Lensing {0}".format(gf))
        gfs = gf.split('/')
        gbase = "/".join(gfs[:-1])
        sfile = "{0}/lensed_{1}".format(sbase, gfs[-1])
        g     = fitsio.read(gf)
        shear = fitsio.read(sfile)

        if 'halo' in gf:
            g     = compute_lensing(g, shear, halos=True)
        else:
            g     = compute_lensing(g, shear)

        gfss = gfs[-1].split('.')
        oname = "{0}/{1}_lensed.{2}.fits".format(gbase, gfss[0], gfss[1])

        fitsio.write(oname, g)


if __name__=='__main__':

  gpath = sys.argv[1]
  sbase = sys.argv[2]

  if '*' in gpath:
      gnames = glob(gpath)
  else:
      gnames = np.loadtxt(gpath)
  
  add_lensing(gnames, sbase)
