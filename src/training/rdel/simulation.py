from __future__ import print_function, division
from helpers import SimulationAnalysis
import scipy.constants as const
import numpy as np
import fitsio

from .abundancematch import abundanceMatchSnapshot

class Simulation(object):

    def __init__(self, boxsize, hfiles, rnnfiles, h, zs=None,
                    zmin=None, zmax=None, zstep=None, nz=None):

        self.boxsize = boxsize
        self.h = h
        self.hfiles = np.array(hfiles)
        self.rnnfiles = np.array(rnnfiles)

        self.lums = np.linspace(-30, -10, 100)

        if zs is not None:
            self.zs = zs
        elif (zmin is not None) & (zmax is not None) & (zstep is not None):
            self.zs = np.arange(zmin, zmax+zstep, zstep)
        elif (zmin is not None) & (zmax is not None) & (nz is not None):
            self.zs = np.linspace(zmin, zmax, nz)

        self.associateFiles()
        self.unitmap = {'mag':'magh', 'phi':'hmpc3dex'}


    def associateFiles(self):

        hfn = np.array([float(h.split('_')[-1].split('.list')[0]) for h in self.hfiles])
        crn = np.array([float(c.split('_')[-1]) for c in self.rnnfiles])

        hidx = hfn.argsort()
        cidx = crn.argsort()

        self.hfiles = self.hfiles[hidx[::-1]]
        self.rnnfiles = self.rnnfiles[cidx[::-1]]

        print(self.hfiles)
        print(self.rnnfiles)
        print(self.zs)

    def getSHAMFileName(self, hfname, alpha, scatter, lfname):

        fs = hfname.split('/')
        fs[-1] = 'sham_{0}_{1}_{2}_{3}'.format(lfname,
                                                alpha,
                                                scatter,
                                                fs[-1])
        fn  = '/'.join(fs)

        return fn

    def abundanceMatch(self, lf, alpha=0.5, scatter=0.17, debug=False,
                       parallel=False, startat=None):
        """
        Abundance match all of the hlists
        """

        odtype = np.dtype([('PX', np.float),
                            ('PY', np.float),
                            ('PZ', np.float),
                            ('AMPROXY', np.float),
                            ('LUMINOSITY', np.float)])

        if parallel:
            from mpi4py import MPI

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            
            hfiles = self.hfiles[rank::size]
            zs     = self.zs[rank::size]
        else:
            hfiles = self.hfiles
            zs = self.zs

        if startat is not None:
            hfiles = hfiles[startat:]

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
            lz = lf.genLuminosityFunctionZ(self.lums, z)
            
            for k in lf.unitmap:
                if lf.unitmap[k] == self.unitmap[k]:
                    continue
                else:
                    conv = self.getattr(self, '{0}2{1}'.format(lf.unitmap[k], self.unitmap[k]))
                    lz[k] = conv(lz[k])

            out['LUMINOSITY'] = abundanceMatchSnapshot(proxy,
                                                        scatter,
                                                        lz,
                                                        self.boxsize,
                                                        debug=debug)

            sfname = self.getSHAMFileName(hf, alpha, scatter,
                                            lf.name)
            fitsio.write(sfname, out)


    def rdelMagDist(self):
        """
        Compute rdel-magnitude distribution in SHAMs
        """

    def mag2magh(mag):
        
        return mag - 5 * np.log10(self.h)

    def mpc3dex2hmpc3dex(phi):

        return phi / self.h ** 3
