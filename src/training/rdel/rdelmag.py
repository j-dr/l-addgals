from __future__ import print_function, division
from scipy.optimize import curve_fit
import cPickle as pickle
import numpy as np

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
    jrmdist, e = np.histogram(rdel, mag, bins=[rbins, magcuts])

    #cumulative brighter than a given magnitude cut
    rmdist     = np.cumsum(jrmdist[:,::-1], axis=1)
    rmerr      = np.sqrt(rmdist)

    #normalize distributions
    rmcounts   = np.sum(rmdist, axis=0).reshape(1,rmdist.shape[1])
    nrmdist    = rmdist / rmcounts
    nrmerr     = rmerr  / rmcounts

    return nrmdist, nrmerr

def fitRdelDist(rdist, rdisterr, rmean):
    """
    Fit a model to a density distribution for a single magnitude cut.
    """

    nmcuts   = rdist.shape[1]
    mphat    = np.zeros(5, nmcuts)
    mpcovhat = np.zeros(5, 5, nmcuts)

    for i in xrange(nmcuts):

        mphat[:,i], mpcovhat[:,:,i] = curve_fit(rdelModel, rmean, rdist,
                                        sigma=rdisterr, p0 = np.ones(5))

    return mphat, mpcovhat

def visRdelSnapshotFit(mphat, mpcovhat, rdist,
                        rdisterr, rmean, mcuts,
                        smname, outdir):

    nmcuts   = len(mcuts)

    xal = int(np.ceil(np.sqrt(nmcuts)))

    f, ax = plt.subplots(xal, xal, sharex=True, sharey=True)

    for i, m in enumerate(mcuts):
        rdisthat = rdelModel(rmean, *mphat[:,i])
        ax[i//xal, i%xal].errorbar(rmean, rdist[:,i], yerr=rdisterr[:,i])
        ax[i//xal, i%xal].plot(rmean, rdisthat)

        chisq = np.sum((rdisthat - rdist[:,i]) ** 2 / rdisterr[:,i]**2)

        ax[i//xal, i%xal].text(1, 1, r"M_{r}<%f.2$" % m)
        ax[i//xal, i%xal].text(1, 0.5, r"$\Chi^{2}=%f.2$" % chisq)

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

    sax.set_xlabel(r'$r_{\delta}\, [Mpc h^{-1}]$', fontsize=18)
    sax.set_ylabel(r'$M_{r}\, [Mag]$', fontsize=18)

    os.makedirs('{}/plots'.format(outdir))

    plt.savefig('{}/plots/{}_rdelmodel.png'.format(outdir, modelname))


def saveRdelSnapshotFit(mphat, mpcovhat, modelname, outdir):

    sdict = {'mphat':mphat, 'mpcovhat':mpcovhat, 'smname':modelname}

    with open('{}/{}_rdelmodel.pkl'.format(outdir, modelname), 'w') as fp:
        pickle.dump(sdict)

def fitSnapshot(shamfile, rnnfile, outdir, debug=False):

    smname = shamfile.split('/')[-1].split('.')
    smname = '.'.join(smname[:-1])

    mag = fitsio.read(shamfile, columns=['LUMINOSITY'])
    rdel = readHaloRdel(rnnfile)

    #set up bins
    rbins = np.logspace(-3., np.log10(15.), 50)
    mcuts = np.linspace(-24, -18, 20)

    rmean = (rbins[1:] + rbins[:-1]) / 2

    rmdist, rmerr = rdelMagDist(rdel, mag)
    mphat, mpcovhat = fitRdelDist(rmdist, rmerr, rmean)

    saveRdelSnapshotFit(mphat, mpcovhat, smname, outdir)

    if debug:
        visRdelSnapshotFit(mphat, mpcovhat, rmdist,
                            rmerr, rmean, mcuts,
                            smname, outdir)

def fitRdelMagZDist(rmzdist):
    pass
