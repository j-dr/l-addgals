from __future__ import print_function, division
from scipy.optimize import curve_fit
from glob import glob
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
    jlcmass  = np.zeros((nmbins, njacktot))

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
            jlcmass[j-1, i] = np.nanmean(lc[midx==j])

    #do the jackknifing

    lcmass    = np.nanmean(jlcmass, axis=1)
    lcmassvar = (np.nanvar(jlcmass, axis=1) *
                    (njacktot - 1))

    return lcmass, lcmassvar, jlcmass


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
            prdisthat = rdelModel(rmean, *pmphat[:,i])
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

    print('Fitting {} and {}'.format(shamfile, rnnfile))

    smname = shamfile.split('/')[-1].split('.')
    smname = '.'.join(smname[:-1])

    mag = fitsio.read(shamfile, columns=['LUMINOSITY'])
    rdel = readHaloRdel(rnnfile)

    #set up bins
    rbins = np.logspace(-3., np.log10(15.), 50)
    mcuts = np.linspace(-24, -18, 20)

    rmean = (rbins[1:] + rbins[:-1]) / 2

    #cludge to temporarily deal with calcrnn leaving out last halo
    if len(mag)==(len(rdel)+1):
        mag = mag[:-1]

    #measure normalized distribution of densities above magnitude cuts
    rmdist, rmerr = rdelMagDist(rdel, mag['LUMINOSITY'], rbins, mcuts)

    #clean up nans and infs
    rmdist, rmerr = cleanDist(rmdist, rmerr)

    #fit models to individual density distributions
    mphat, mpcovhat = fitRdelDist(rmdist, rmerr, rmean)

    saveRdelSnapshotFit(mphat, mpcovhat, rmdist, rmerr, smname, outdir)

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
    scales = []
    for m in modelfiles:
        ms = m.split('hlist_')[-1]
        mss = ms.split('.list')
        scales.append(float(mss[0]))
        
    scales = np.array(scales)
    sidx   = scales.argsort()
    z      = 1 / scales[sidx[::-1]] - 1
    models = models[sidx[::-1]]

    return models, z

def loadModels(inbase, smbase):

    models = []
    szs    = []
    zs     = []
    for i in range(len(smbase)):
        smodels, szm = loadSnapshotFits(inbase, smbase[i])
        models.append(smodels)
        szs.append(set(szm))
        zs.append(szm)

    z = np.array(list(set.intersection(*szs)))
    z.sort()

    for i in range(len(models)):
        zi = np.array(zs[i])
        i1d = np.in1d(zi, z)
        models[i] = models[i][i1d]
        zidx = zi[i1d].argsort()
        models[i] = models[i][zidx]

    return models, z


def fitRdelMagZDist(inbase, smbase, z, mcuts, zmcut,
                    zmin=0.0, zmax=2.0, fmlim=-18, 
                    bmlim=-22.5, validate_fit=True):

    #load single snapshot models
    models, z = loadModels(inbase, smbase)

    #construct parameter arrays
    #select parameter fits from appropriate simulation for each z-mag bin
    parr, pearr, mucarr, mucearr, sigmacarr, sigmacearr, mufarr, mufearr, sigmafarr, sigmafearr = makeOutputArrays(models, zmcut, z, mcuts)

    mag_ref = -20.5
    xvec, zmidx, zmaidx, mbidx, mfidx = makeVandermonde(z, mcuts, bmlim, 
                                                        fmlim, mag_ref, 
                                                        zmin, zmax)

    #construct x, y vectors to predict for
    xpvec, _, _, _, _ = makeVandermonde(z, mcuts, np.min(mcuts)-1, np.max(mcuts)+1, 
                            mag_ref, np.min(z)-1, np.max(z)+1)

    pf, mucf, muff, sigmacf, sigmaff = setupOutputs(parr, mucarr, sigmacarr, mufarr,
                                                    sigmafarr, zmidx, zmaidx, mbidx, mfidx)

    betap, betamuc, betamuf, betasigmacf, betasigmaff = solveLeastSquares(xvec, pf,
                                                                          mucf, muff,
                                                                          sigmacf, sigmaff)

    phatarr, muchatarr, sigmachatarr, mufhatarr, sigmafhatarr = makePrediction(xpvec, betap,
                                                                               betamuc,
                                                                               betamuf,
                                                                               betasigmacf,
                                                                               betasigmaff,
                                                                               len(z), len(mcuts))

    rbins = np.logspace(-3., np.log10(15.), 50)
    rmean = (rbins[1:] + rbins[:-1]) / 2

    if validate_fit:
        validateRdelMagZDist(phatarr, muchatarr, sigmachatarr,
                             mufhatarr, sigmafhatarr, models,
                             rmean, mcuts, z, inbase)

    saveRdelMagZDistFit(inbase, smbase, betap[0], betamuc[0], betamuf[0],
                        betasigmacf[0], betasigmaff[0])

    return betap, betamuc, betamuf, betasigmacf, betasigmaff, models

def makeOutputArrays(models, zmcut, z, mcuts):

    parr = np.zeros((len(z), len(models[0][0]['mphat'][0,:])))
    pearr = np.zeros((len(z), len(models[0][0]['mphat'][0,:])))
    mucarr = np.zeros((len(z), len(models[0][0]['mphat'][1,:])))
    mucearr = np.zeros((len(z), len(models[0][0]['mphat'][1,:])))
    sigmacarr = np.zeros((len(z), len(models[0][0]['mphat'][2,:])))
    sigmacearr = np.zeros((len(z), len(models[0][0]['mphat'][2,:])))
    mufarr = np.zeros((len(z), len(models[0][0]['mphat'][3,:])))
    mufearr = np.zeros((len(z), len(models[0][0]['mphat'][3,:])))
    sigmafarr = np.zeros((len(z), len(models[0][0]['mphat'][4,:])))
    sigmafearr = np.zeros((len(z), len(models[0][0]['mphat'][4,:])))

    for i in xrange(len(z)):
        for j in xrange(len(mcuts)):
            midx = zmcut[i,j]
            parr[i,j] = models[midx][i]['mphat'][0,j]
            pearr[i,j] = models[midx][i]['mpcovhat'][0,0,j]
            mucarr[i,j] = models[midx][i]['mphat'][1,j]
            mucearr[i,j] = models[midx][i]['mpcovhat'][1,1,j]
            sigmacarr[i,j] = models[midx][i]['mphat'][2,j]
            sigmacearr[i,j] = models[midx][i]['mpcovhat'][2,2,j]
            mufarr[i,j] = models[midx][i]['mphat'][3,j]
            mufearr[i,j] = models[midx][i]['mpcovhat'][3,3,j]
            sigmafarr[i,j] = models[midx][i]['mphat'][4,j]
            sigmafearr[i,j] = models[midx][i]['mpcovhat'][4,4,j]

    return parr, pearr, mucarr, mucearr, sigmacarr, sigmacearr, mufarr, mufearr, sigmafarr, sigmafearr


def makePrediction(x, betap, betamuc, betamuf, betasigmac, betasigmaf, nz, nm):
    phat = np.dot(x.T, betap[0])
    muchat = np.dot(x.T, betamuc[0])
    mufhat = np.dot(x.T, betamuf[0])
    sigmachat = np.dot(x.T, betasigmac[0])
    sigmafhat = np.dot(x.T, betasigmaf[0])

    phatarr = np.array(phat).reshape(nz, nm)
    muchatarr = np.array(muchat).reshape(nz, nm)
    mufhatarr = np.array(mufhat).reshape(nz, nm)
    sigmachatarr = np.array(sigmachat).reshape(nz, nm)
    sigmafhatarr = np.array(sigmafhat).reshape(nz, nm)

    return phatarr, muchatarr, mufhatarr, sigmachatarr, sigmafhatarr

def solveLeastSquares(x, p, muc, muf, sigmac, sigmaf):

    betap = np.linalg.lstsq(x.T, p)
    betamuc = np.linalg.lstsq(x.T, muc)
    betamuf = np.linalg.lstsq(x.T, muf)
    betasigmacf = np.linalg.lstsq(x.T, sigmac)
    betasigmaff = np.linalg.lstsq(x.T, sigmaf)

    return betap, betamuc, betamuf, betasigmacf, betasigmaff
    
def makeVandermonde(z, mcuts, bmlim, fmlim, mag_ref, zmin, zmax):

    mag_ref = -20.5
    bright_mag_lim = bmlim-mag_ref
    faint_mag_lim  = fmlim-mag_ref

    zmidx  = z.searchsorted(zmin)
    zmaidx = z.searchsorted(zmax)

    mbidx = mcuts.searchsorted(bmlim)
    mfidx = mcuts.searchsorted(fmlim)

    #construct x, y vectors to fit to
    x = np.meshgrid(mcuts[mbidx:mfidx], z[zmidx:zmaidx])
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

    return xvec, zmidx, zmaidx, mbidx, mfidx

def setupOutputs(parr, mucarr, sigmacarr, mufarr, sigmafarr, zmidx, zmaidx, mbidx, mfidx):

    pf   = parr[zmidx:zmaidx, mbidx:mfidx].flatten().reshape(-1, 1)
    mucf = mucarr[zmidx:zmaidx, mbidx:mfidx].flatten().reshape(-1,1)
    sigmacf = sigmacarr[zmidx:zmaidx, mbidx:mfidx].flatten().reshape(-1,1)
    muff = mufarr[zmidx:zmaidx, mbidx:mfidx].flatten().reshape(-1,1)
    sigmaff = sigmafarr[zmidx:zmaidx, mbidx:mfidx].flatten().reshape(-1,1)

    return pf, mucf, sigmacf, muff, sigmaff


def getCoeffNames(rank):
    
    if rank==4:
        names = ['z0', 'z1', 'z2', 'z3', 'z4', 'm1',
                 'm2', 'm3', 'm4', 'm1z1', 'm1z2',
                 'm2z1', 'mz2', 'm1z3', 'm2z2',
                 'm3z1']

    return names


def saveRdelMagZDistFit(outbase, smbase, betap, betamuc, betamuf, 
                          betasigmac, betasigmaf, rank=4):

    names = getCoeffNames(rank)

    with open('{}/{}_global_rdelmodel.txt'.format(outbase, smbase), 'w') as fp:
        for i, name in enumerate(names):
            fp.write('p{} {}\n'.format(name, betap[i]))
        for i, name in enumerate(names):
            fp.write('cm{} {}\n'.format(name, betamuc[i]))
        for i, name in enumerate(names):
            fp.write('cs{} {}\n'.format(name, betasigmac[i]))
        for i, name in enumerate(names):
            fp.write('fm{} {}\n'.format(name, betamuf[i]))
        for i, name in enumerate(names):
            fp.write('fs{} {}\n'.format(name, betasigmaf[i]))



def validateRdelMagZDist(parr, mcarr, scarr, mfarr,
                          sfarr, inmodels, rmean, 
                          mcuts, z, outdir):

    pmphat = np.zeros((5, len(mcuts), len(z)))
    for i, zi in enumerate(z):
        for j, mc in enumerate(mcuts):
            pmphat[:,j,i] = [parr[i,j], mcarr[i,j], scarr[i,j],
                                 mfarr[i,j], sfarr[i,j]]

    for i, zi in enumerate(z):
        f, ax = visRdelSnapshotFit(inmodels[0][i]['mphat'],
                            inmodels[0][i]['mpcovhat'],
                            inmodels[0][i]['rmdist'],
                            inmodels[0][i]['rmerr'],
                            rmean, mcuts,
                            'pred_'+inmodels[0][i]['smname'],
                            outdir, pmphat=pmphat[:,:,i])
