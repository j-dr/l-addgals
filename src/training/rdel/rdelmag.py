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

def lcenModel(mass, m0, mc, a, b, k):

    return (mr0 - 2.5 * ( a * np.log10(mass / mc) -
              np.log10(1 + (mass/mc) ** (b * k)) / k))


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

def lcenMassDist(x, y, z, lcen, mass, mbins, lbox, njackside=5):

    nmbins   = len(mbins) - 1
    njacktot = njackside ** 3
    jlcmass  = np.zeros((nmbins, njackside))

    #want jackknife errors on lcenmass
    xi = (njackside * x) // lbox
    yi = (njackside * y) // lbox
    zi = (njackside * z) // lbox

    bidx = xi * njackside ** 2 + yi * njackside + zi

    for i in xrange(njacktot):
        idx = bidx == i
        lc  = lcen[idx]
        m   = mass[idx]

        midx = np.digitize(m, mbins)

        for j in xrange(1, nmbins+1):
            jlcmass[j, i] = np.mean(lc[midx==j])

    #do the jackknifing

    lcmass    = np.mean(jlcmass, axis=1)
    lcmassvar = (np.sum(jlcmass - lcmass.reshape(nmbins, 1), axis=1) *
                    (njacktot - 1) / njacktot)

    return lcmass, lcmassvar


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
                        smname, outdir, pmphat=None):

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

        if pmphat is not None:
            prdisthat = rdelmodel(rmean, *pmphat[:,i])
            ax[i//xal, i%xal].plot(rmean, prdisthat)

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

def saveRdelSnapshotFit(mphat, mpcovhat, rmdist, rmerr, modelname, outdir):

    sdict = {'mphat':mphat, 'mpcovhat':mpcovhat, 'rmdist':rmdist, 'rmerr':rmerr, 'smname':modelname}

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

    saveRdelSnapshotFit(mphat, mpcovhat, rmdist, smname, outdir)

    if debug:
        visRdelSnapshotFit(mphat, mpcovhat, rmdist,
                            rmerr, rmean, mcuts,
                            smname, outdir)

def loadSnapshotFits(inbase, smbase):

    modelfiles = glob('{}/{}*rdelmodel.pkl'.format(inbase, smbase))
    models     = []

    for m in modelfiles:
        with open(m, 'r') as fp:
            models.append(pickle.load(fp))

    models = np.array(models)
    scales = np.array([float(m.split('_')[-2]) for m in modelfiles])
    sidx   = scales.argsort
    z      = 1 / scales[sidx[::-1]] - 1
    models = models[sidx[::-1]]

    return models, z

def fitRdelMagZDist(inbase, smbase, z, mcuts, zmcut,
                zmin=0.0, zmax=2.0, fmlim=-18, bmlim=-22.5):

    #load single snapshot models
    models = []

    for i in range(len(smbase)):
        smodels, szm = loadSnapshotFits(inbase, smbase[i])
        models.append(smodels)



    #construct parameter arrays
    parr = np.zeros((len(z), len(models[0]['mphat'][0,:])))
    pearr = np.zeros((len(z), len(models[0]['mphat'][0,:])))
    mucarr = np.zeros((len(z), len(models[0]['mphat'][1,:])))
    mucearr = np.zeros((len(z), len(models[0]['mphat'][1,:])))
    sigmacarr = np.zeros((len(z), len(models[0]['mphat'][2,:])))
    sigmacearr = np.zeros((len(z), len(models[0]['mphat'][2,:])))
    mufarr = np.zeros((len(z), len(models[0]['mphat'][3,:])))
    mufearr = np.zeros((len(z), len(models[0]['mphat'][3,:])))
    sigmafarr = np.zeros((len(z), len(models[0]['mphat'][4,:])))
    sigmafearr = np.zeros((len(z), len(models[0]['mphat'][4,:])))

    #select parameter fits from appropriate simulation for each z-mag bin
    for i in xrange(len(z)):
        for j in xrange(len(mcuts)):
            midx = zmcut[i,j]
            parr[i,j] = models[midx]['mphat'][0,j]
            pearr[i,j] = models[midx]['mpcovhat'][0,0,j]
            mucarr[i,j] = models[midx]['mphat'][1,j]
            mucearr[i,j] = models[midx]['mpcovhat'][1,1,j]
            sigmacarr[i,j] = models[midx]['mphat'][2,j]
            sigmacearr[i,j] = models[midx]['mpcovhat'][2,2,j]
            mufarr[i,j] = models[midx]['mphat'][3,j]
            mufearr[i,j] = models[midx]['mpcovhat'][3,3,j]
            sigmafarr[i,j] = models[midx]['mphat'][4,j]
            sigmafearr[i,j] = models[midx]['mpcovhat'][4,4,j]

    mag_ref = -20.5
    bright_mag_lim = bmlim-mag_ref
    faint_mag_lim  = fmlim-mag_ref

    zmidx  = z.searchsorted(zmin)
    zmaidx = z.searchsorted(zmax)

    mbidx = mcuts.searchsorted(bmlim)
    mfidx = mcuts.searchsorted(fmlim)

    #construct x, y vectors to fit to
    x = np.meshgrid(mcuts[bmidx:mfidx], z[zmidx:zmaidx])
    zv = x[1].flatten() -1.0
    mv = x[0].flatten()
    mv = mv - mag_ref
    mv[mv<bright_mag_lim]  = bright_mag_lim
    mv[mv>faint_mag_lim]   = faint_mag_lim
    o  = np.ones(len(zv))

    #construct vandermonde matrix
    xvec = np.array([o, zv, zv**2, zv**3, zv**4,
                        mv, mv**2, mv**3, mv**4,
                        mv*zv, mv*zv**2, mv**2*zv,
                        (mv*zv)**2, mv*zv**3, mv**2*zv**2,
                         mv**3*zv])

    pf   = parr[zmidx:zmaidx, bmidx:bfidx].flatten().reshape(-1, 1)
    mucf = mucarr[zmidx:zmaidx, bmidx:bfidx].flatten().reshape(-1,1)
    sigmacf = sigmacarr[zmidx:zmaidx, bmidx:bfidx].flatten().reshape(-1,1)
    muff = mufarr[zmidx:zmaidx, bmidx:bfidx].flatten().reshape(-1,1)
    sigmaff = sigmafarr[zmidx:zmaidx, bmidx:bfidx].flatten().reshape(-1,1)

    #construct x, y vectors to predict for
    x = np.meshgrid(mcuts, z)
    zv = x[1].flatten() -1.0
    mv = x[0].flatten()
    mv = mv - mag_ref
    mv[mv<bright_mag_lim]  = bright_mag_lim
    mv[mv>faint_mag_lim]   = faint_mag_lim
    o  = np.ones(len(zv))

    xpvec = np.array([o, zv, zv**2, zv**3,
                        zv**4, mv, mv**2, mv**3,
                        mv**4, mv*zv, mv*zv**2,
                        mv**2*zv, (mv*zv)**2,
                        mv*zv**3, mv**2*zv**2,
                        mv**3*zv])

    betap = np.linalg.lstsq(xvec.T, pf)
    betamuc = np.linalg.lstsq(xvec.T, mucf)
    betamuf = np.linalg.lstsq(xvec.T, muff)
    betasigmacf = np.linalg.lstsq(xvec.T, sigmacf)
    betasigmaff = np.linalg.lstsq(xvec.T, sigmaff)
    phat = np.dot(xpvec.T, betap[0])
    muchat = np.dot(xpvec.T, betamuc[0])
    mufhat = np.dot(xpvec.T, betamuf[0])
    sigmachat = np.dot(xpvec.T, betasigmacf[0])
    sigmafhat = np.dot(xpvec.T, betasigmaff[0])
    phatarr = np.array(phat).reshape(zidx, len(mcuts)-1)
    muchatarr = np.array(muchat).reshape(zidx, len(mcuts)-1)
    mufhatarr = np.array(mufhat).reshape(zidx, len(mcuts)-1)
    sigmachatarr = np.array(sigmachat).reshape(zidx, len(mcuts)-1)
    sigmafhatarr = np.array(sigmafhat).reshape(zidx, len(mcuts)-1)

    validateRdelMagZDist(phatarr, muchatarr, sigmachatarr,
                            mufhatarr, sigmafhararr, models,
                            mcuts, z, outdir)


def validateRdelMagZDist(parr, mcarr, scarr, mfarr,
                          sfarr, inmodels, rmdist,
                          rmean, rmdisterr, mcuts,
                          z, outdir):

    pmphat = np.zeros((5, len(mcuts), len(z)))
    for i, zi in enumerate(z):
        for j, mc in enumerate(mcuts):
            pmphat[:,j,i] = [parr[i,j], mcarr[i,j], scarr[i,j],
                                 mfarr[i,j], sfarr[i,j]]

    for i, zi in enumerate(z):
        f, ax = visRdelSnapshotFit(inmodels[i]['mphat'],
                            inmodels[i]['mpcovhat'],
                            inmodels[i]['rmdist'],
                            inmodels[i]['rmerr'],
                            rmean, mcuts,
                            'pred_'+models[i]['smname'],
                            outdir, pmphat=pmphat)
