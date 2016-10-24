from __future__ import print_function, division
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
import numpy as np
import fitsio
import os

from .io import readHaloRdel


def rdelModel(rmean, p, muc, sigmac, muf, sigmaf):
    """
    Evaluate a gaussian + lognormal distribution
    """

    return (1 - p) * np.exp(-(np.log(rmean) - muc) ** 2 / (2 * sigmac ** 2)) / ( rmean * np.sqrt(2 * np.pi ) * sigmac ) + p * np.exp(-(rmean - muf) ** 2 / (2 * sigmaf ** 2)) / (np.sqrt(2 * np.pi ) * sigmaf )

def rdelMagDist(rdel, mag, rbins, magcuts):
    """
    Calculate the distribution of rdel of galaxies brighter
    than some absolute magniutude cut
    """

    #calculate the joint distribution
    jrmdist, e0, e1 = np.histogram2d(rdel, mag, bins=[rbins, magcuts])

    #cumulative brighter than a given magnitude cut
    rmdist     = np.cumsum(jrmdist, axis=1)
    rmerr      = np.sqrt(rmdist)

    #normalize distributions
    rmcounts   = np.sum(rmdist, axis=0).reshape(1,rmdist.shape[1])
    nrmdist    = rmdist / rmcounts / (rbins[1:] - rbins[:-1]).reshape(rbins.shape[0]-1,1)
    nrmerr     = rmerr  / rmcounts / (rbins[1:] - rbins[:-1]).reshape(rbins.shape[0]-1,1)

    return nrmdist, nrmerr

def cleanDist(rmdist, rmerr):

    rmdist[np.isnan(rmdist) | ~np.isfinite(rmdist)] = 0.0
    rmerr[np.isnan(rmdist) | ~np.isfinite(rmdist)]  = 100000.0
    rmerr[rmerr==0.] = 0.001

    return rmdist, rmerr

def fitRdelDist(rdist, rdisterr, rmean, useerror=False):
    """
    Fit a model to a density distribution for a single magnitude cut.
    """

    nmcuts   = rdist.shape[1]
    mphat    = np.zeros((5, nmcuts))
    mpcovhat = np.zeros((5, 5, nmcuts))

    for i in xrange(nmcuts):
        if i==0:
            p0 = np.array([0.661101, -0.654770, 0.554604, 2.58930, 1.00034])
        else:
            p0 = mphat[:,nmcuts-i]

        try:

            if useerror:
                mphat[:,nmcuts-i-1], mpcovhat[:,:,nmcuts-i-1] = curve_fit(rdelModel, 
                                                                          rmean, 
                                                                          rdist[:,nmcuts-i-1],
                                                                          sigma=rdisterr[:,nmcuts-i-1], 
                                                                          p0 = p0, 
                                                                          absolute_error=True)
            else:
                mphat[:,nmcuts-i-1], mpcovhat[:,:,nmcuts-i-1] = curve_fit(rdelModel, 
                                                                          rmean, 
                                                                          rdist[:,nmcuts-i-1],
                                                                          p0 = p0)
        except RuntimeError as e:
            print('Fit failed for magnitude cut: {0}'.format(nmcuts-i))
            mphat[:,nmcuts-i-1] = p0
            mpcovhat[:,:,nmcuts-i-1] = mpcovhat[:,:,nmcuts-i-1]

    return mphat, mpcovhat

def visRdelSnapshotFit(mphat, mpcovhat, rdist,
                        rdisterr, rmean, mcuts,
                        smname, outdir):

    nmcuts   = len(mcuts)

    xal = int(np.ceil(np.sqrt(nmcuts)))

    f, ax = plt.subplots(xal, xal, sharex=True, sharey=False)

    sax = f.add_subplot(111)
    plt.setp(sax.get_xticklines(), visible=False)
    plt.setp(sax.get_yticklines(), visible=False)
    plt.setp(sax.get_xticklabels(), visible=False)
    plt.setp(sax.get_yticklabels(), visible=False)
    sax.patch.set_alpha(0.0)
    sax.patch.set_facecolor('none')
    sax.spines['top'].set_color('none')
    sax.spines['bottom'].set_color('none')
    sax.spines['left'].set_color('none')
    sax.spines['right'].set_color('none')
    sax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    sax.set_xlabel(r'$R_{\delta}\, [Mpc h^{-1}]$', fontsize=18, labelpad=20)
    sax.set_ylabel(r'$M_{r}\, [Mag]$', fontsize=18, labelpad=20)

    for i, m in enumerate(mcuts[:-1]):
        rdisthat = rdelModel(rmean, *mphat[:,i])
        ax[i//xal, i%xal].errorbar(rmean, rdist[:,i], yerr=rdisterr[:,i])
        ax[i//xal, i%xal].plot(rmean, rdisthat)

        ax[i//xal, i%xal].set_xscale('log')
        
        chisq = np.sum((rdisthat - rdist[:,i]) ** 2 / rdisterr[:,i]**2)

        #ax[i//xal, i%xal].text(0.1, 0.9, r"$M_{r}<%.2f$" % m, fontsize=14)
                    
        #ax[i//xal, i%xal].text(0.1, 0.7, r"$\chi^{2}=%.2e$" % chisq, fontsize=14)

    try:
        os.makedirs('{}/plots'.format(outdir))
    except:
        pass

    f.set_figheight(15)
    f.set_figwidth(15)

    plt.savefig('{}/plots/{}_rdelmodel.png'.format(outdir, smname))

    return f, ax

def saveRdelSnapshotFit(mphat, mpcovhat, modelname, outdir):

    sdict = {'mphat':mphat, 'mpcovhat':mpcovhat, 'smname':modelname}

    with open('{}/{}_rdelmodel.pkl'.format(outdir, modelname), 'w') as fp:
        pickle.dump(sdict, fp)

def fitSnapshot(shamfile, rnnfile, outdir, debug=False):

    smname = shamfile.split('/')[-1].split('.')
    smname = '.'.join(smname[:-1])

    mag = fitsio.read(shamfile, columns=['LUMINOSITY'])
    rdel = readHaloRdel(rnnfile)

    #set up bins
    rbins = np.logspace(-3., np.log10(15.), 50)
    mcuts = np.linspace(-24, -18, 20)

    rmean = (rbins[1:] + rbins[:-1]) / 2

    #measure normalized distribution of densities above magnitude cuts
    rmdist, rmerr = rdelMagDist(rdel, mag['LUMINOSITY'], rbins, mcuts)
    
    #clean up nans and infs
    rmdist, rmerr = cleanDist(rmdist, rmerr)
    
    #fit models to individual density distributions
    mphat, mpcovhat = fitRdelDist(rmdist, rmerr, rmean)

    saveRdelSnapshotFit(mphat, mpcovhat, smname, outdir)

    if debug:
        visRdelSnapshotFit(mphat, mpcovhat, rmdist,
                            rmerr, rmean, mcuts,
                            smname, outdir)

def fitRdelMagZDist(rmzdist):
    pass
