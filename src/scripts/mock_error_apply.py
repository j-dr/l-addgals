from __future__ import print_function, division
from glob import glob
from mpi4py import MPI
import numpy.lib.recfunctions as rf
import numpy as np
import healpy as hp
import fitsio
import pickle
import yaml
import sys
import os

models = {
    'DR8':
        {'maglims':[20.425,21.749,21.239,20.769,19.344],
        'exptimes' : [21.00,159.00,126.00,99.00,15.00],
        'lnscat' : [0.284,0.241,0.229,0.251,0.264]
         },

    'STRIPE82':
        {
        'maglims' : [22.070,23.432,23.095,22.649,21.160],
        'exptimes' : [99.00,1172.00,1028.00,665.00,138.00],
        'lnscat' : [0.280,0.229,0.202,0.204,0.246]
        },

    'CFHTLS':
        {
        'maglims' : [24.298,24.667,24.010,23.702,22.568],
        'exptimes' : [2866.00,7003.00,4108.00,3777.00,885.00],
        'lnscat' : [0.259,0.244,0.282,0.258,0.273]
        },
    'DEEP2':
        {
        'maglims' : [24.730,24.623,24.092],
        'exptimes' : [7677.00,8979.00,4402.00],
        'lnscat' : [0.300,0.293,0.300]
        },
    'FLAMEX':
        {
        'maglims' : [21.234,20.929],
        'exptimes' : [259.00,135.00],
        'lnscat' : [0.300,0.289]
        },
    'IRAC':
        {
        'maglims' : [19.352,18.574],
        'exptimes' : [8.54,3.46],
        'lnscat' : [0.214,0.283]
        },
    'NDWFS':
        {
        'maglims' : [25.142,23.761,23.650],
        'exptimes' : [6391.00,1985.00,1617.00],
        'lnscat' : [0.140,0.294,0.272]
        },

    'RCS':
        {
        'maglims' : [23.939,23.826,23.067,21.889],
        'exptimes' : [2850.00,2568.00,1277.00,431.00],
        'lnscat' : [0.164,0.222,0.250,0.271]
        },
    'VHS':
        {
        'maglims' : [20.141,19.732,19.478],
        'exptimes' : [36.00,31.00,23.00],
        'lnscat' : [0.097,0.059,0.069]
        },
    'VIKING':
        {
        'maglims' : [21.643,20.915,20.768,20.243,20.227],
        'exptimes' : [622.00,246.00,383.00,238.00,213.00],
        'lnscat' : [0.034,0.048,0.052,0.040,0.066]
        },
    'DC6B':
        {
        'maglims' : [24.486,23.473,22.761,22.402],
        'exptimes' : [2379.00,1169.00,806.00,639.00],
        'lnscat' : [0.300,0.300,0.300,0.300]
        },

    'DES':
        {
        'maglims' : [24.956,24.453,23.751,23.249,21.459],
        'exptimes' : [14467.00,12471.00,6296.00,5362.00,728.00],
        'lnscat' : [0.2,0.2,0.2,0.2,0.2]
        },

    'BCS_LO':
        {
        'maglims' : [23.082,22.618,22.500,21.065],
        'exptimes' : [809.00,844.00,641.00,108.00],
        'lnscat' : [0.277,0.284,0.284,0.300]
        },

    'BCS':
        {
        'maglims' : [23.360,23.117,22.539,21.335],
        'exptimes' : [838.00,1252.00,772.00,98.00],
        'lnscat' : [0.276,0.272,0.278,0.279]
        },
    
    'DES_SV': 
    {
	'maglims' : [23.621,23.232,23.008,22.374,20.663],
	'exptimes' : [4389.00,1329.00,1405.00,517.00,460.00],
	'lnscat' : [0.276,0.257,0.247,0.241,0.300] 
        },
    
    'DES_SV_OPTIMISTIC':
        {
        'maglims' : [23.621+0.5,23.232+0.5,23.008,22.374,20.663],
        'exptimes' : [4389.00,1329.00,1405.00,517.00,460.00],
        'lnscat' : [0.276,0.257,0.247,0.241,0.300]
        },
    'WISE':
        {
        'maglims' : [19.352,18.574],
        'exptimes' : [8.54,3.46],
        'lnscat' : [0.214,0.283]
        },

    'DECALS':
        {
       'maglims' : [23.3,23.3,22.2,20.6,19.9],
       'exptimes' : [1000,3000,2000,1500,1500],
       'lnscat' : [0.2,0.2,0.2,0.2,0.2]
       },

}

def calc_nonuniform_errors(exptimes,limmags,mag_in,nonoise=False,zp=22.5,nsig=10.0,
                           fluxmode=False,lnscat=None,b=None,inlup=False,detonly=False):

    seed = time()

    f1lim = 10**((limmags-zp)/(-2.5))
    fsky1 = (((f1lim**2)*exptimes)/(nsig**2.) - f1lim)
    fsky1[fsky1<0.001] = 0.001

    if inlup:
        bnmgy = b*1e9
        tflux = exptimes*2.0*bnmgy*np.sinh(-np.log(b)-0.4*np.log(10.0)*mag_in)
    else:
        tflux = exptimes*10**((mag_in - zp)/(-2.5))
        
    noise = np.sqrt(fsky1*exptimes + tflux)

    if nonoise:
        flux = tflux
    else:
        flux = tflux + noise*np.random.randn(len(mag_in))

    if lnscat is not None:
        noise = np.exp(np.log(noise) + lnscat*np.random.randn(len(mag_in)))

    if fluxmode:
        mag = flux/exptimes
        mag_err = noise/exptimes
    else:
        if b is not None:
            bnmgy = b*1e9
            flux_new = flux/exptimes
            noise_new = noise/exptimes
            mag = 2.5*np.log10(1.0/b) - asinh2(0.5*flux_new/(bnmgy))/(0.4*np.log(10.0))

            mag_err = 2.5*noise_new/(2.*bnmgy*np.log(10.0)*np.sqrt(1.0+(0.5*flux_new/(bnmgy))**2.))

        else:
            mag = zp-2.5*np.log10(flux/exptimes)
            mag_err = (2.5/np.log(10.))*(noise/flux)

            bad = np.where(np.isfinite(mag)==False)
            nbad = len(bad)
            if (detonly):
                mag[bad]=99.0
                mag_err[bad]=99.0
                
    return mag, mag_err


def calc_uniform_errors(model, tmag, maglims, exptimes, lnscat):


    nmag=len(maglims)
    ngal=len(tmag)

    zp=22.5
    #calculate fsky1 -- sky in 1 second
    flux1_lim = 10**((maglims-zp)/(-2.5)) 
    flux1_lim[flux1_lim < 120/exptimes] = 120/exptimes[flux1_lim < 120/exptimes]
    fsky1 = (flux1_lim**2*exptimes)/100. - flux1_lim
    
    oflux=np.zeros((ngal, nmag))
    ofluxerr=np.zeros((ngal, nmag))
    omag=np.zeros((ngal, nmag))
    omagerr=np.zeros((ngal, nmag))
    offset = 0.0

    for i in range(nmag):
        tflux = exptimes[i] * 10**((tmag[:,i]-offset-zp)/(-2.5))
        noise = np.exp(np.log(np.sqrt(fsky1[i]*exptimes[i] + tflux)) 
                    + lnscat[i]*np.random.randn(ngal))

        flux = tflux + noise*np.random.randn(ngal)
        
        oflux[:,i] = flux / exptimes[i]
        ofluxerr[:,i] = noise/exptimes[i]
        
        omag[:,i] = 22.5-2.5*np.log10(oflux[:,i])
        omagerr[:,i] = (2.5/np.log(10.))*(ofluxerr[:,i]/oflux[:,i])

        bad,=np.where(~np.isfinite(omag[:,i]))
        nbad = len(bad)
        if (nbad > 0) :
            omag[bad,i] = 99.0
            omagerr[bad,i] = 99.0

    
    return omag, omagerr, oflux, ofluxerr

def apply_nonuniform_errormodel(fname, obase, depthfile, magfile=None, usemag=None,
                                survey=None, rot=None, rotcols=None):
    
    g = fitsio.read(fname)

    if rot is not None:
        vec = np.dot(rot, g[['PX','PY','PZ']].view((g['PX'].dtype,3)).T)
        
        if rotcols is not None:
            adtype = [np.dtype([(f,np.float)]) for f in rotcols]
            theta, phi = hp.vec2ang(vec[:,0], vec[:,1], vec[:,2])
                                                     
            dec = 90-theta*180./np.pi
            ra  = phi*180./np.pi
            
            if 'DEC' in rotcols[0]:
                data = [dec, ra]
            else:
                data = [ra, dec]

                g = rf.append_fields(g,rotcols, data=data,
                                     dtypes=adtype, usemask=False)

    if magfile is not None:
        mags = fitsio.read(magfile)
        if ('LMAG' in mags.dtype.names) and (mags['LMAG']!=0).any():
            imtag = 'LMAG'
            omag = mags['LMAG']
        else:
            imtag = 'TMAG'
            omag = mags['TMAG']
        fs = magfile.split('/')
        oname = "{0}/{1}".format(obase,magfile)
    else:
        if ('LMAG' in mags.dtype.names) and (mags['LMAG']!=0).any():
            imtag = 'LMAG'
            omag = g['LMAG']
        else:
            imtag = 'TMAG'
            omag = g['TMAG']

        fs = fname.split('/')
        fss = fname.split('.')
        oname = "{0}/{1}_{2}.{3}.fits".format(obase, fss[:-2], survey, fss[-2])

    #get mags to use
    if usemag is None:
        nmag = omag.shape[1]
        usemag = range(nmag)
    else:
        nmag = len(usemag)

    #make output structure
    if magfile is not None:
        obs = mags
    else:
        obs = np.zeros(len(g), dtype=np.dtype([(imtag,(np.float,nmag)),('AMAG',(np.float,nmag)),
                                               ('OMAG',(np.float,nmag)),('OMAGERR',(np.float,nmag)),
                                               ('FLUX',(np.float,nmag)),('IVAR',(np.float,nmag))]))
        for i,idx in enumerate(usemag):
            obs[imtag][:,i] = g[imtag][:,idx]
        
    d, dhdr = fitsio.read(depthfile, header=True)

    if (survey=="Y1A1") | (survey=="DES") | (survey=="SVA"):
        mindec = -90.
        maxdec = 90
        minra = 0.0
        maxra = 360.

    elif survey=="DR8":
        mindec = -20
        maxdec = 90
        minra = 0.0
        maxra = 360.

    maxtheta=(90.0-mindec)*np.pi/180.
    mintheta=(90.0-maxdec)*np.pi/180.
    minphi=minra*np.pi/180.
    maxphi=maxra*np.pi/180.
    
    #keep pixels in footprint
    theta, phi = hp.pix2ang(dhdr['NSIDE'],d['HPIX'])
    infp = np.where(((mintheta < theta) & (theta < maxtheta)) & ((minphi < phi) & (phi < maxphi)))
    d = d[infp]

    #match galaxies to correct pixels of depthmap
    pixind = d['HPIX'].argsort()
    d = d[pixind]

    if rot is None:
        theta = (90-g['DEC'])*np.pi/180.
        phi   = (g['RA']*np.pi/180.)
        pix   = hp.ang2pix(dhdr['NSIDE'],theta, phi)
    else:
        pix   = hp.ang2pix(dhdr['NSIDE'],vec[:,0], vec[:,1], vec[:,2])

    pixind = d['HPIX'].searchsorted(pix)
    guse = (ipring!=0)&(ipring!=len(d))

    if not any(guse):
        print("No galaxies in this pixel are in the footprint")
        return

    for ind,i in enumerate(usemag):
        
        flux, fluxerr = calc_nonuniform_errors(d['EXPTIMES'][pixind[guse],ind],
                                               d['LIMMAGS'][pixind[guse],ind],
                                               mag_in[guse,i], fluxmode=True)
        obs['FLUX'][guse,ind] = flux
        obs['IVAR'][guse,ind] = 1/fluxerr**2
        obs['OMAG'][guse,ind] = 22.5 - 2.5*np.log10(flux)
        obs['OMAGERR'][guse,ind] = 1.086*fluxerr/flux

        obs['FLUX'][~guse,ind] = -99
        obs['IVAR'][~guse,ind] = -99
        obs['OMAG'][~guse,ind] = 99
        obs['OMAGERR'][~guse,ind] = 99

        bad = (flux<=0)

        obs['OMAG'][guse[bad],ind] = 99.0
        obs['OMAGERR'][guse[bad],ind] = 99.0


        r = np.random.rand(len(pixind[guse]))
        bad, = np.where(r>depth['FRACGOOD'][pixind[guse]])

        if len(bad)>0:
            obs['OMAG'][guse[bad],ind] = 99.0
            obs['OMAGERR'][guse[bad],ind] = 99.0

    
    fitsio.write(oname, obs)

    if rotcols is not None:
        fitsio.write(fname+'.test', g)

        
    

def apply_uniform_errormodel(fname, obase, model, detonly=False, magname=None, usemags=None):

    print("Working on {0}".format(fname))
    g = fitsio.read(fname)

    maglims = np.array(models[model]['maglims'])
    exptimes = np.array(models[model]['exptimes'])
    lnscat = np.array(models[model]['lnscat'])

    omag, omagerr, oflux, ofluxerr = calc_uniform_errors(model, g['TMAG'],
                                                         maglims, exptimes,
                                                         lnscat)
    print('Done calculating errors')
    nmags = omag.shape[1]

    emdtype = np.dtype([('OMAG', (omag.dtype, nmags)), 
                        ('OMAGERR',(omagerr.dtype,nmags)),
                        ('OFLUX', (oflux.dtype, nmags)),
                        ('OFLUXERR', (ofluxerr.dtype, nmags))])

    fs = fname.split('.')
    fp = fs[-2]
    pb = fs[-4]
    pb = pb.split('/')[-1]

    oname = "{0}/{1}_{2}.{3}.fits".format(obase, pb, model, fp)
    print('Oname: {0}'.format(oname))

    out = np.zeros(len(omag), dtype=emdtype)
    out['OMAG']     = omag
    out['OMAGERR']  = omagerr
    out['OFLUX']    = oflux
    out['OFLUXERR'] = ofluxerr

    if detonly:
        idx = omag[omag<99.]
        idx = idx.any(axis=1)
        out = out[idx]

    fitsio.write(oname, out)

if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    cfgfile = sys.argv[1]
    with open(cfgfile, 'r') as fp:
        cfg = yaml.load(fp)

    gpath = cfg['GalPath']
    model = cfg['Model']
    obase = cfg['OutBase']

    fnames = np.array(glob(gpath))
    
    if 'DepthFile' in cfg.keys():
        dfile = cfg['DepthFile']
        uniform = False
    else:
        uniform = True

    if ('MagPath' in cfg.keys()) and (cfg['MagPath'] is not None):
        mnames  = np.array(glob(cfg['MagPath']))

        fpix = np.array([int(f.split('.')[-2]) for f in fnames])
        mpix = np.array([int(f.split('.')[-2]) for f in mnames])

        fidx = fpix.argsort()
        midx = mpix.argsort()

        assert((fpix[fidx]==mpix[midx]).all())

        fnames = fnames[fidx]
        mnames = mnames[midx]

    else:
        mnames = [None]*len(fnames)

    if 'Rotation' in cfg.keys():
        rfile = cfg['Rotation']
        with open(rfile, 'r') as fp:
            rot   = pickle.load(fp)
        if 'RotCols' in cfg.keys():
            rotcols = cfg['RotCols']
        else:
            rotcols = None
            

    else:
        rot = None

    if 'UseMags' in cfg.keys():
        usemags = cfg['UseMags']

    if rank==0:
        try:
            os.makedirs(obase)
        except Exception as e:
            print(e)

    for fname, mname in zip(fnames[rank::size],mnames[rank::size]):
        if uniform:
            apply_uniform_errormodel(fname, obase, model, magfile=mname, usemags=usemags)
        else:
            apply_nonuniform_errormodel(fname,obase,dfile,magfile=mname,survey=model,rot=rot, rotcols=rotcols)
