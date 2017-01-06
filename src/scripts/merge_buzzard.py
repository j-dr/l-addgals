import numpy as np
import fitsio as fio
import healpy as hp
import glob
import yaml
import os
import sys

'''
'''

class buzzard_flat_cat(object):

    def __init__(
        self,
        rootdir   = '/project/projectdirs/des/jderose/addgals/catalogs/Buzzard/Catalog_v1.1/', 
        obsdir    = 'des_obs_rotated/Y1A1/', 
        pzdir     = 'des_obs_rotated/Y1A1/bpz/',
        truthdir  = 'truth_rotated/',
        obsname   = 'Y1A1.',
        truthname = 'truth.',
        pzname    = 'Y1A1_bpz.',
        debug     = False,
        simnum    = 0):

        self.maxrows=500000000

        self.rootdir   = rootdir
        self.truthdir  = truthdir
        self.truthname = truthname
        self.obsdir    = obsdir
        self.obsname   = obsname
        self.pzdir     = pzdir
        self.pzname    = pzname
        self.simnum    = str(simnum)

        self.loop_cats(debug=debug)

        return

    def loop_cats(self, debug=False): 

        if debug:
            gold  = np.zeros(self.maxrows, dtype = [('coadd_objects_id','i8')] 
                             + [('ra','f4')]
                             + [('dec','f4')]
                             + [('redshift','f4')]
                             + [('mag_g', 'f4')]
                             + [('mag_r', 'f4')]
                             + [('mag_i', 'f4')]
                             + [('mag_z', 'f4')]
                             + [('mag_y', 'f4')]
                             + [('flags_badregion','i8')] 
                             + [('flags_gold','i8')]
                             + [('hpix','i8')]
                             + [('sample','i8')])
        else:
            gold  = np.zeros(self.maxrows, dtype = [('coadd_objects_id','i8')] 
                             + [('ra','f4')]
                             + [('dec','f4')]
                             + [('redshift','f4')]
                             + [('flags_badregion','i8')] 
                             + [('flags_gold','i8')]
                             + [('hpix','i8')]
                             + [('sample','i8')])            

        if debug:
            shape = np.zeros(self.maxrows, dtype = [('coadd_objects_id','i8')] 
                             + [('e1','f4')] 
                             + [('e2','f4')] 
                             + [('g1','f4')] 
                             + [('g2','f4')] 
                             + [('kappa','f4')] 
                             + [('m1','f4')] 
                             + [('m2','f4')] 
                             + [('c1','f4')] 
                             + [('c2','f4')] 
                             + [('weight','f4')]
                             + [('size','f4')]
                             + [('flags','i8')])
        else:
            shape = np.zeros(self.maxrows, dtype = [('coadd_objects_id','i8')] 
                             + [('e1','f4')] 
                             + [('e2','f4')] 
                             + [('g1','f4')] 
                             + [('g2','f4')] 
                             + [('kappa','f4')] 
                             + [('m1','f4')] 
                             + [('m2','f4')] 
                             + [('c1','f4')] 
                             + [('c2','f4')] 
                             + [('weight','f4')]
                             + [('flags','i8')])


        photoz = np.zeros(self.maxrows, dtype = [('coadd_objects_id', 'i8')]
                                          + [('mean_z', 'f8')]
                                          + [('z_mc', 'f8')]    
                                          + [('redshift', 'f8')]    
                                          + [('weight', 'f8')]    
                                          + [('flags', 'f8')])

        lenst = 0

        gout = fio.FITS('Buzzard_v1.1_'+self.simnum+'_gold.fits.gz', 'rw') 
        sout = fio.FITS('Buzzard_v1.1_'+self.simnum+'_shape.fits.gz', 'rw') 
        pout = fio.FITS('Buzzard_v1.1_'+self.simnum+'_pz.fits.gz', 'rw')

        for ifile,filename in enumerate(glob.glob(self.rootdir+self.obsdir+'*'+self.obsname+'*.fit')):
            tname = filename.replace(self.obsname, self.truthname).replace(self.obsdir, self.truthdir)
            pzname = filename.replace(self.obsname, self.pzname).replace(self.obsdir, self.pzdir).replace('fit', 'fits')

            truth  = fio.FITS(tname)[-1].read(columns=['ID','GAMMA1','GAMMA2','KAPPA','Z'])
            if not debug:
                obs    = fio.FITS(filename)[-1].read(columns=['RA','DEC','EPSILON1','EPSILON2','LSS_FLAG', 'WL_FLAG'])
            else:
                obs    = fio.FITS(filename)[-1].read(columns=['RA','DEC','EPSILON1','EPSILON2','LSS_FLAG', 'WL_FLAG', 'MAG_G', 'MAG_R', 'MAG_I', 'MAG_Z', 'MAG_Y', 'SIZE'])
            pz     = fio.FITS(pzname)[-1].read(columns=['MEAN_Z','Z_MC'])

            sflag = obs['LSS_FLAG'] + 2 * obs['WL_FLAG']
            truth = truth[sflag>0]
            obs   = obs[sflag>0]
            pz    = pz[sflag>0]
            sflag = sflag[sflag>0]

            # insert selection function here to mask truth/obs (if can be run on individual files)

            print ifile, len(truth), filename

            gold['coadd_objects_id'][lenst:lenst+len(truth)]  = truth['ID']
            gold['ra'][lenst:lenst+len(truth)]                = obs['RA']
            gold['dec'][lenst:lenst+len(truth)]               = obs['DEC']
            gold['redshift'][lenst:lenst+len(truth)]          = truth['Z']
            gold['hpix'][lenst:lenst+len(truth)]              = hp.ang2pix(4096, np.pi/2.-np.radians(obs['DEC']),np.radians(obs['RA']), nest=True)
            gold['sample'][lenst:lenst+len(truth)]            = sflag
            if debug:
                gold['mag_g'][lenst:lenst+len(truth)]         = obs['MAG_G']
                gold['mag_r'][lenst:lenst+len(truth)]         = obs['MAG_R']
                gold['mag_i'][lenst:lenst+len(truth)]         = obs['MAG_I']
                gold['mag_z'][lenst:lenst+len(truth)]         = obs['MAG_Z']
                gold['mag_y'][lenst:lenst+len(truth)]         = obs['MAG_Y']

            shape['coadd_objects_id'][lenst:lenst+len(truth)] = truth['ID']
            shape['e1'][lenst:lenst+len(truth)]               = obs['EPSILON1']
            shape['e2'][lenst:lenst+len(truth)]               = obs['EPSILON2']
            shape['g1'][lenst:lenst+len(truth)]               = truth['GAMMA1']
            shape['g2'][lenst:lenst+len(truth)]               = truth['GAMMA2']
            shape['kappa'][lenst:lenst+len(truth)]            = truth['KAPPA']
            shape['m1'][lenst:lenst+len(truth)]               += 1.
            shape['m2'][lenst:lenst+len(truth)]               += 1.
            shape['weight'][lenst:lenst+len(truth)]           += 1.
            if debug:
                shape['size'][lenst:lenst+len(truth)]         = obs['SIZE']

            photoz['coadd_objects_id'][lenst:lenst+len(truth)] = truth['ID']
            photoz['mean_z'][lenst:lenst+len(truth)]           = pz['MEAN_Z']
            photoz['z_mc'][lenst:lenst+len(truth)]             = pz['Z_MC']
            photoz['redshift'][lenst:lenst+len(truth)]         = truth['Z']
            photoz['weight'][lenst:lenst+len(truth)]           += 1.

            if ifile == 0:
                gout.write(gold[lenst:lenst+len(truth)])
                sout.write(shape[lenst:lenst+len(truth)])
                pout.write(photoz[lenst:lenst+len(truth)])
            else:
                gout[-1].append(gold[lenst:lenst+len(truth)])
                sout[-1].append(shape[lenst:lenst+len(truth)])
                pout[-1].append(photoz[lenst:lenst+len(truth)])

            lenst+=len(truth)
        
        gout.close()
        sout.close()
        pout.close()

        return 

if __name__ == '__main__':
    
    cfgfile = sys.argv[1]

    with open(cfgfile, 'r') as fp:
        cfg = yaml.load(fp)

    obj = buzzard_flat_cat(**cfg)
