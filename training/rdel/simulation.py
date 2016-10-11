from __future__ import print_function, division
from helpers import SimulationAnalysis
import numpy as np
import fitsio

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


    def getSHAMFileName(hfname, alpha, scatter, lfname):

        fs = hfname.split('/')
        fs[-1] = 'sham_{0}_{1}_{2}_{3}'.format(lfname,
                                                alpha,
                                                scatter,
                                                fs[-1])
        fn  = '/'.join(fs)

        return fn

    def abundanceMatch(lf, alpha=0.5, scatter=0.17):
        """
        Abundance match all of the hlists
        """

        lfz      = np.zeros(len(self.lums), 2)
        lfz[:,0] = self.lums

        odtype = np.dtype([('PX', np.float),
                            ('PY', np.float),
                            ('PY', np.float),
                            ('AMPROXY', np.float),
                            ('LUMINOSITY', np.float)])

        for i, hf in enumerate(self.hfiles):

            halos = SimulationAnalysis.readHlist(hf,
                                                    ['vmax',
                                                     'vvir',
                                                     'x',
                                                     'y',
                                                     'z'])
            proxy = halos['vvir'] * (halos['vmax'] / halos['vvir']) ** alpha

            out = np.zeros(len(proxy), dtype=odtype)

            out['PX'] = halos['x']
            out['PY'] = halos['y']
            out['PZ'] = halos['z']
            out['AMPROXY'] = proxy

            z = self.zs[i]
            lfz[:,1] = lf.genLuminosityFunctionZ(self.lums, z)

            out['LUMINOSITY'] = abundanceMatchSnapshot(proxy,
                                                        scatter,
                                                        lfz,
                                                        self.boxsize)

            sfname = self.getSHAMFileName(hf, alpha, scatter,
                                            lf.name)
            fitsio.write(sfname, out)


    def rdelMagDist():
        """
        Compute rdel-magnitude distribution in SHAMs
        """
