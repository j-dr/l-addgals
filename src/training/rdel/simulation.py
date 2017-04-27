from __future__ import print_function, division
from helpers import SimulationAnalysis
from itertools import izip
import scipy.constants as const
import numpy as np
import warnings
import fitsio
import os

from .abundancematch import abundanceMatchSnapshot
from .rdelmag        import fitSnapshot

class Simulation(object):

    def __init__(self, name, boxsize, snapdirs, 
                   hfiles, rnnfiles, outdir, 
                   h, zs=None, zmin=None, zmax=None, 
                   zstep=None, nz=None, shamfiles=None,
                   compressed_hlist=False, simtype='LGADGET2',
                   scaleh=False, scaler=False, strscale=None):

        self.name     = name
        self.boxsize  = boxsize
        self.h        = h
        self.snapdirs= np.array(snapdirs)
        self.hfiles   = np.array(hfiles)
        self.rnnfiles = np.array(rnnfiles)
        self.outdir   = outdir
        self.simtype  = simtype
        self.compressed_hlist = compressed_hlist
        self.scaleh = scaleh
        self.scaler = scaler
        self.strscale  = strscale

        if shamfiles is not None:
            self.shamfiles = np.array(shamfiles)
        else:
            self.shamfiles = shamfiles

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

        hfn = np.array([int(h.split('_')[-1].split('.')[0]) for h in self.hfiles])
        hz = self.zs[hfn]
        shz = self.strscale[hfn]

        crn = np.array([int(c.split('_')[-1]) for c in self.rnnfiles])
        cz = self.zs[crn]


        if len(hfn) > len(crn):
            inz = np.in1d(hfn, crn)
            self.hfiles = self.hfiles[inz]
            hfn = hfn[inz]
            hz  = hz[inz]
            shz = shz[inz]
            inz = np.in1d(crn, hfn)
            self.rnnfiles = self.rnnfiles[inz]
            crn = crn[inz]
            cz  = cz[inz]

        else:
            inz = np.in1d(crn, hfn)
            self.rnnfiles = self.rnnfiles[inz]
            crn = crn[inz]
            cz  = cz[inz]
            inz = np.in1d(hfn, crn)
            self.hfiles = self.hfiles[inz]
            hfn = hfn[inz]
            hz  = hz[inz]
            shz = shz[inz]

        assert(len(self.hfiles)==len(self.rnnfiles))

        hidx = hfn.argsort()
        cidx = crn.argsort()

        self.hfiles   = self.hfiles[hidx[::-1]]
        self.rnnfiles = self.rnnfiles[cidx[::-1]]
        self.zs       = hz[hidx[::-1]]
        self.strscale    = shz[hidx[::-1]]
        self.nums     = hfn[hidx[::-1]]

        print(self.hfiles)
        print(self.rnnfiles)
        print(self.zs)

    def getSHAMFileName(self, hfname, alpha, scatter, lfname, ind):

        fs = hfname.split('/')

        if self.nums[ind]<10:
            num  = '00{}'.format(self.nums[ind])
        else:
            num = '0{}'.format(self.nums[ind])

        fs[-1] = 'sham_{}_{}_{}_{}_{}'.format(self.name, 
                                                lfname,
                                                alpha,
                                                scatter,
                                                fs[-1])

        print(fs[-1])
        print(num)
        fn = fs[-1].replace(num, '{}.list'.format(self.strscale[ind]))

        return fn

    def abundanceMatch(self, lf, alpha=0.7, scatter=0.17, debug=False,
                       parallel=False, startat=None):
        """
        Abundance match all of the hlists
        """

        odtype = np.dtype([('PX', np.float),
                            ('PY', np.float),
                            ('PZ', np.float),
                            ('AMPROXY', np.float),
                            ('VMAX', np.float),
                            ('MVIR', np.float),
                            ('RVIR', np.float),
                            ('LUMINOSITY', np.float),
                            ('CENTRAL', np.int)])

        if parallel:
            from mpi4py import MPI

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()

            hfiles = self.hfiles[rank::size]
            zs     = self.zs[rank::size]

        else: 
            rank = 0
            hfiles = self.hfiles
            zs = self.zs

        if startat is not None:
            hfiles = hfiles[startat:]
            zs = zs[startat:]
            
        if rank == 0:
            try:
                os.makedirs('{0}/sham/'.format(self.outdir))
            except OSError as e:
                warnings.warn('Directory {0}/sham/ already exists!'.format(self.outdir), Warning)
                pass

        for i, hf in enumerate(hfiles):
            if startat is not None:
                ind = startat + rank + i * size
            else:
                ind = startat + rank + i * size

            sfname = self.getSHAMFileName(hf, alpha, scatter,
                                                  lf.name, ind)

            oname = '{0}/sham/{1}'.format(self.outdir, sfname)

            if os.path.isfile(oname):
                print('{} exists! Skipping this snapshot'.format(oname))
                continue 

            if not self.compressed_hlist:
                halos = SimulationAnalysis.readHlist(hf,
                                                     ['vmax',
                                                      'mvir',
                                                      'rvir',
                                                      'upid',
                                                      'x',
                                                      'y',
                                                      'z'])
            else:
                halos = fitsio.read(hf, columns=['vmax', 'mvir',
                                                 'rvir',
                                                 'upid',
                                                 'x','y','z'])

            vvir  = (const.G * halos['mvir'] / halos['rvir']) ** 0.5
            proxy = vvir * (halos['vmax'] / vvir) ** alpha

            out = np.zeros(len(proxy), dtype=odtype)

            out['PX'] = halos['x']
            out['PY'] = halos['y']
            out['PZ'] = halos['z']
            out['VMAX'] = halos['vmax']
            out['MVIR'] = halos['mvir']
            out['RVIR'] = halos['rvir']
            out['CENTRAL'][halos['upid']==-1] = 1
            out['CENTRAL'][halos['upid']!=-1] = 0
            out['AMPROXY'] = proxy

            z = zs[i]
            lz = lf.genLuminosityFunctionZ(self.lums, z)

            for k in lf.unitmap:
                if lf.unitmap[k] == self.unitmap[k]:
                    continue
                else:
                    conv = getattr(self, '{0}2{1}'.format(lf.unitmap[k], self.unitmap[k]))
                    lz[k] = conv(lz[k])

            out['LUMINOSITY'] = abundanceMatchSnapshot(proxy,
                                                        scatter,
                                                        lz,
                                                        self.boxsize,
                                                        debug=debug)

            try:
                fitsio.write('{0}/sham/{1}'.format(self.outdir, sfname), out)
            except IOError as e:
                print('File {} already exists, not writing new one!'.format('{0}/sham/{1}'.format(self.outdir, sfname)))
                

    def getSHAMFiles(self, lf, alpha=0.7, scatter=0.17):
        
        shamfiles = []
        
        for i, hf in enumerate(self.hfiles):
            shamfiles.append('{}/sham/{}'.format(self.outdir, self.getSHAMFileName(hf, alpha, scatter, lf.name, i)))

        self.shamfiles = np.array(shamfiles)

    def rdelMagDist(self, lf, debug=False, 
                      startat=None, parallel=False,
                      alpha=0.7, scatter=0.17):
        """
        Compute rdel-magnitude distribution in SHAMs
        """

        if self.shamfiles is None:
            self.getSHAMFiles(lf, alpha=alpha, scatter=scatter)

        if parallel:
            from mpi4py import MPI

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()

            shamfiles = self.shamfiles[rank::size]
            rnnfiles  = self.rnnfiles[rank::size]
        else:
            shamfiles = self.shamfiles
            rnnfiles  = self.rnnfiles
            rank = 0

        if startat is not None:
            shamfiles = shamfiles[startat:]
            rnnfiles = rnnfiles[startat:]

        if rank == 0:
            try:
                os.makedirs('{}/rdel'.format(self.outdir))
            except OSError as e:
                warnings.warn('Directory {}/rdel already exists!', Warning)

        for sf, rf in izip(shamfiles, rnnfiles):
            fitSnapshot(sf, rf, '{}/rdel/'.format(self.outdir), debug=debug)

        #read in all models, fit global model


    def mag2magh(self, mag):

        return mag - 5 * np.log10(self.h)

    def mpc3dex2hmpc3dex(self, phi):

        return phi / self.h ** 3
