from __future__ import print_function, division
from helpers import SimulationAnalysis
import scipy.constants as const
import numpy as np
import fitsio

if __name__ == '__main__':
    import mpi4py as MPI

from .abundancematch import abundanceMatchSnapshot

class Simulation(object):

    def __init__(self, boxsize, hfiles, rnnfiles, zs=None,
                    zmin=None, zmax=None, zstep=None, nz=None):

        self.boxsize = boxsize
        self.hfiles = hfiles
        self.rnnfiles = rnnfiles

        self.lums = np.linspace(-27, -10, 100)

        if zs is not None:
            self.zs = zs
        elif (zmin is not None) & (zmax is not None) & (zstep is not None):
            self.zs = np.arange(zmin, zmax+zstep, zstep)
        elif (zmin is not None) & (zmax is not None) & (nz is not None):
            self.zs = np.linspace(zmin, zmax, nz)


    def associateFiles(self):

        hfn = np.array([h.split('_')[-1].split('.list')[0] for h in self.hfiles])
        crn = np.array([c.split('_')[-1] for c in self.rnnfiles])

        hidx = hfn.argsort()
        cidx = crn.argsort()

        self.hfiles = self.hfiles[hidx]
        self.rnnfiles = self.rnnfiles[cidx]

    def getSHAMFileName(self, hfname, alpha, scatter, lfname):

        fs = hfname.split('/')
        fs[-1] = 'sham_{0}_{1}_{2}_{3}'.format(lfname,
                                                alpha,
                                                scatter,
                                                fs[-1])
        fn  = '/'.join(fs)

        return fn

    def abundanceMatch(self, lf, alpha=0.5, scatter=0.17, debug=False,
                       parallel=False):
        """
        Abundance match all of the hlists
        """

        lfz      = np.zeros((len(self.lums), 2))
        lfz[:,0] = self.lums

        odtype = np.dtype([('PX', np.float),
                            ('PY', np.float),
                            ('PZ', np.float),
                            ('AMPROXY', np.float),
                            ('LUMINOSITY', np.float)])

        if parallel:
            comm = comm.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            
            hfiles = self.hfiles[rank::size]
            zs     = self.zs[rank::size]
        else:
            hfiles = self.hfiles
            

        for i, hf in enumerate(hfiles):

            halos = SimulationAnalysis.readHlist(hf,
                                                    ['vmax',
                                                     'mvir',
                                                     'rvir',
                                                     'x',
                                                     'y',
                                                     'z'])

            vvir  = (const.G * halos['mvir'] / halos['rvir']) ** 0.5
            proxy = vvir * (halos['vmax'] / vvir) ** alpha

            out = np.zeros(len(proxy), dtype=odtype)

            out['PX'] = halos['x']
            out['PY'] = halos['y']
            out['PZ'] = halos['z']
            out['AMPROXY'] = proxy

            z = zs[i]
            lfz[:,1] = lf.genLuminosityFunctionZ(self.lums, z)

            out['LUMINOSITY'] = abundanceMatchSnapshot(proxy,
                                                        scatter,
                                                        lfz,
                                                        self.boxsize,
                                                        debug=debug)

            sfname = self.getSHAMFileName(hf, alpha, scatter,
                                            lf.name)
            fitsio.write(sfname, out)


    def rdelMagDist(self):
        """
        Compute rdel-magnitude distribution in SHAMs
        """
