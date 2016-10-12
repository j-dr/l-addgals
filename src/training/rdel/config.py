from __future__ import print_function, division
from glob import glob
import numpy as np
import yaml

from .simulation import Simulation
from .model      import Model
from .luminosityfunction import LuminosityFunction, DSGLuminosityFunction


def readCfg(filename):

    with open(filename, 'r') as fp:
        cfg = yaml.load(fp)

    return cfg

def parseConfig(cfg):

    simcfg = cfg['Simulation']
    sims = []

    for i, s in enumerate(simcfg['hlistbase']):
        hlists = glob('{0}/hlist*list'.format(s))
        rnn    = glob('{0}/snapdir_*/rnn*[0-9]'.format(simcfg['rnnbase']))
        if 'snaptimes' in simcfg:
            a = np.loadtxt(simcfg['snaptimes'][i])
            zs = 1/a[:,1] - 1.
            sims.append(Simulation(simcfg['boxsize'][i],
                                   hlists,
                                   rnn,
                                   zs=zs))

        else:
            sims.append(Simulation(simcfg['boxsize'][i],
                                   hlists,
                                   rnn,
                                   zmin=simcfg['zmin'][i],
                                   zmax=simcfg['zmax'][i],
                                   nz=simcfg['nz'][i]))

    lf = DSGLuminosityFunction()
    m  = Model(sims)

    return sims, lf, m
