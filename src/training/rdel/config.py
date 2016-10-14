from __future__ import print_function, division
from glob import glob
import numpy as np
import yaml

from .simulation import Simulation
from .model      import Model
from .luminosityfunction import LuminosityFunction, DSGLuminosityFunction, BernardiLuminosityFunction


def readCfg(filename):

    with open(filename, 'r') as fp:
        cfg = yaml.load(fp)

    return cfg

def setLF(cfg):
    
    if cfg['LuminosityFunction']['type'] == 'DSG':

        lf = DSGLuminosityFunction()
    
    elif cfg['LuminosityFunction']['type'] == 'Bernardi':
        
        lf = BernardiLuminosityFunction(cfg['LuminosityFunction']['Q'])

    return lf
                       
        

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
                                   simcfg['h'][i],
                                   zs=zs))

        else:
            sims.append(Simulation(simcfg['boxsize'][i],
                                   hlists,
                                   rnn,
                                   simcfg['h'][i],
                                   zmin=simcfg['zmin'][i],
                                   zmax=simcfg['zmax'][i],
                                   nz=simcfg['nz'][i]))

    lf = setLF(cfg)
    m  = Model(sims)

    return sims, lf, m
