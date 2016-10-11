from __future__ import print_function, division
from glob import glob

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
        hlists = glob('{0}/hlist*'.format(s))
        rnn    = glob('{0}/rnn*'.format(simcfg['rnnbase']))

        sims.append(Simulation(simcfg['boxsize'][i],
                                hlists,
                                rnn,
                                zmin=simcfg['zmin'][i],
                                zmax=simcfg['zmax'][i],
                                nz=simcfg['nz'])

    lf = DSGLuminosityFunction()
    m  = Model(sims)
    
