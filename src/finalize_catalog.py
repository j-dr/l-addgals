from __future__ import print_function, division
from SimulationAnalysis import readHlist
from collections import deque
from glob import glob
from mpi4py import MPI
from scipy import spatial
import numpy.lib.recfunctions as rf
import numpy as np
import healpy as hp
import fitsio
import time


TZERO = None
def tprint(info):
    global TZERO
    if TZERO is None:
        TZERO = time.time()

    print('[%8ds] %s' % (time.time()-TZERO,info))


def finalize_catalogs(basepath, prefix, postfix, outpath, halopaths, mmin=[3e12,3e12,2.4e13]):

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    tags = {'write':0, 'occ1050':1, 'occ2600':2, 'occ4000':3,
            'lum1050':4, 'lum2600':5, 'lum4000':6, 'exit':7}
    message = None

    if rank == 0:
        writing = deque()
        wwaiting = deque()
        done = []
        occ = [None, None, None]
        lum = [None, None, None] 
        while True:
            #are we done?
            if len(done)==(size-1): break
            remove = []
            #see if any write requests can be filled
            if len(wwaiting)>0:

                for w in wwaiting:
                    if w[1] not in writing:
                        comm.send(message, tag=tags['write'], dest=w[0])
                        writing.append(w[1])
                        remove.append(w)
            
            for r in remove:
                wwaiting.remove(r)

            #Recieve requests
            status = MPI.Status()
            message = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            tag = status.Get_tag()
            if tag==tags['write']:
                wwaiting.append([status.Get_source(), message])
            elif tag==tags['occ1050']:
                if occ[0]==None:
                    occ[0] = message
                else:
                    occ[0] += message
            elif tag==tags['occ2600']:
                if occ[1]==None:
                    occ[1] = message
                else:
                    occ[1] += message
            elif tag==tags['occ4000']:
                if occ[2]==None:
                    occ[2] = message
                else:
                    occ[2] += message
            elif tag==tags['lum1050']:
                if lum[0]==None:
                    lum[0] = message
                else:
                    lum[0] += message
            elif tag==tags['lum2600']:
                if lum[1]==None:
                    lum[1] = message
                else:
                    lum[1] += message
            elif tag==tags['lum4000']:
                if lum[2]==None:
                    lum[2] = message
                else:
                    lum[2] += message
            elif tag==tags['exit']:
                done.append(status.Get_source())
        
        for i, bsize in enumerate(['1050', '2600', '4000']):
            update_halo_file(halopaths[i], prefix, outpath, bsize, occ[i], lum[i], mmin[i])
            
    if rank != 0:
        for i, bsize in enumerate(['1050', '2600', '4000']):
            tprint('    {0}: Processing box {1}'.format(rank, bsize))
            h = fitsio.read(halopaths[i], columns=['HALOID', 'MVIR', 'RVIR', 'HALOPX', 'HALOPY', 'HALOPZ', 'HOST_HALOID'])
            h = h[h['MVIR']>mmin[i]]
            occ = np.zeros((len(h),6))
            lum = np.zeros((len(h),3))
            
            pixh = hp.vec2pix(2, h['HALOPX'], h['HALOPY'], h['HALOPZ'])

            #index the halos by healpix cell
            pidx = pixh.argsort()
            pixh = pixh[pidx]
            h = h[pidx]
            pidx = pixh[1:]-pixh[:-1]
            pidx, = np.where(pidx!=0)
            pidx = np.hstack([np.zeros(1,dtype=np.int),pidx+1])
            upix = pixh[pidx]
            pixmap = hp.ud_grade(np.arange(12*2**2),4)

            pixpaths = glob('{0}/Lb{1}_{2}/[0-9]*'.format(basepath, bsize, postfix))
            chunks = [pixpaths[i::(size-1)] for i in range(size-1)]

            for ppath in chunks[rank-1]:
                #build kdtree for this cell
                pix = int(ppath.split('/')[-1])
                hidx, = np.where(upix==pixmap[pix])
                hstart = pidx[hidx]
                if hidx==(len(pidx)-1):
                    hend = len(h)
                else:
                    hend = pidx[hidx+1]
                
                tprint('    {0}: Building KDTree for pixel {1}'.format(rank,pix))
                ht = spatial.KDTree(h[hstart:hend][['HALOPX','HALOPY','HALOPZ']].view((h.dtype['HALOPX'],3)))
                pfiles = glob('{0}/*/hv_output/gal_ginfo1.dat'.format(ppath))
                
                for i, f in enumerate(pfiles):
                    pc, hdr = fitsio.read(f, header=True)
                    o, l = associate_halos(pc, h[hstart:hend], ht)
                    occ[hstart:hend] += o
                    lum[hstart:hend] += l
                    pc['ID'] += i*1000000000

                    #request permission to write
                    status = MPI.Status()
                    comm.send(pix, tag=tags['write'])
                    message = comm.recv(tag=tags['write'], status=status)
                    write_fits(outpath, prefix, pix, pc, hdr)
            
            comm.send(occ, tag=tags['occ'+bsize])
            comm.send(lum, tag=tags['lum'+bsize])

        comm.send(message, tag=tags['exit'])

def associate_halos(galaxies, halos, tree):

    occ = np.zeros((len(halos),6))
    lum = np.zeros((len(halos),3))
    
    #assign nearest halo to galaxy
    d, hid = tree.query(galaxies[['PX', 'PY', 'PZ']].view((galaxies.dtype['PX'],3)))
    print('Distance to nearest halo: {0}'.format(d))
    print('Minimum distance to nearest halo: {0}'.format(np.min(d)))

    cen = galaxies['CENTRAL']==1
    galaxies['RHALO'][cen] = 0.0
    galaxies['RHALO'][~cen] = 999999
    galaxies['RHALO'][~cen] = d[~cen]
    galaxies['HALOID'][~cen] = halos['HALOID'][hid[~cen]]
    galaxies['M200'][~cen] = halos['MVIR'][hid[~cen]]
    bound = galaxies['RHALO']<=halos['RVIR'][hid]
    lum[hid[bound],0] += galaxies['AMAG'][bound,1]

    nmt, = np.where(halos['HALOID'][hid[cen]]==galaxies[cen]['HALOID'])
    print(len(nmt))
    print('Number of centrals with unidentified halos: {0}'.format(len(galaxies[cen])-len(nmt)))
    galaxies[cen]['M200'] = halos[hid[cen]]['MVIR']
    lum[hid[cen],2] = galaxies[cen]['AMAG'][:,1]

    occ[hid[bound],0] += 1

    for i, mr in enumerate([-10,-18,-19,-20,-21,-22]):
        midx = galaxies['AMAG'][bound,1]<mr
        hidx = hid[bound][midx].argsort()
        mhid = hid[bound][midx][hidx]
        mgr = galaxies['AMAG'][bound,1][midx][hidx]
        hidx = mhid[1:]-mhid[:-1]
        hidx, = np.where(hidx!=0)
        hidx = np.hstack([np.zeros(1,dtype=np.int),hidx+1])
        uhid = mhid[hidx]
        hidx = np.hstack([hidx, np.array([len(mhid)])])
        if (mr==-10):
            for j, uid in enumerate(uhid):
                lum[uid,0]+=np.sum(mgr[hidx[j]:hidx[j+1]])
            continue
        elif (mr==-20):
            for j, uid in enumerate(uhid):
                lum[uid,1]+=np.sum(mgr[hidx[j]:hidx[j+1]])
        
        occ[uhid,i] += hidx[1:]-hidx[:-1]

    return occ, lum

def write_fits(outpath, prefix, pix, catalog, header, nside=8):

    hpix = hp.vec2pix(nside, catalog['PX'],catalog['PY'],catalog['PZ'])
    upix = np.unique(hpix)

    for p in upix:
        fits = fitsio.FITS('{0}/{1}_{2}.fits'.format(outpath, prefix, p) ,'rw')
        try:
            pidx = hpix==p
            h = fits[-1].read_header()
            fits[-1].write_key('NAXIS2', h['NAXIS2']+len(catalog[pidx]))
            fits[-1].append(catalog[pidx])
        except:
            fits.write(catalog, header=header)

def update_halo_file(halopath, prefix, outpath, bsize, occ, lum, mmin):
    
    h = fitsio.read(halopath)
    h = h[h['MVIR']>mmin]
    h['LUMTOT'] = lum[:,0]
    h['LUM20'] = lum[:,1]
    h['LBCG'] = lum[:,2]
    h['NGALS'] = occ[:,0]
    h['N18'] = occ[:,1]
    h['N19'] = occ[:,2]
    h['N20'] = occ[:,3]
    h['N21'] = occ[:,4]
    h['N22'] = occ[:,5]

    fitsio.write('{0}/{1}_halos_{2}.fits'.format(outpath,prefix,bsize), h)

def combine_parents_list_fits(parentpath, listpath, prefix, outpath, bsize, occ, lum, mmin=5e12):
    hdtype = np.dtype([('id',np.int), ('mvir',np.float), ('vmax',np.float), ('vrms',np.float),
                       ('rvir',np.float), ('rs',np.float), ('np',np.int), ('x',np.float), ('y',np.float),
                       ('z',np.float), ('vx',np.float), ('vy',np.float), ('vz',np.float), ('jx',np.float),
                       ('jy',np.float), ('jz',np.float), ('spin',np.float), ('rs_klypin',np.float),
                       ('mvir_all',np.float), ('m200b',np.float), ('m200c',np.float), ('m500c',np.float),
                       ('m2500c',np.float), ('xoff',np.float), ('voff',np.float), ('lambda',np.float),
                       ('b_to_a',np.float), ('c_to_a',np.float), ('ax',np.float), ('ay',np.float),
                       ('az',np.float), ('virial_ratio', np.float)])
    usecols = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
    print('Reading parents')
    parents = np.genfromtxt(parentpath, usecols=[2,14])
    parents = parents[parents[:,0]>mmin]
    print('Reading hlist')
    hlist = np.loadtxt(listpath, usecols=usecols, dtype=hdtype)
    hlist = hlist[hlist['mvir']>mmin]
    print('Making data')
    adtype = np.dtype([('lumtot',np.float), ('lum20',np.float), ('lbcg', np.float),
                       ('ngals',np.int), ('n18',np.int), ('n19',np.int), ('n20',np.int),
                       ('n21',np.int), ('n22',np.int), ('host_haloid',np.int)])
    data = np.zeros(len(hlist), dtype=adtype)
    data['lumtot'] = lum[:,0]
    data['lum20'] = lum[:,1]
    data['lbcg'] = lum[:,2]
    data['ngals'] = occ[:,0]
    data['n18'] = occ[:,1]
    data['n19'] = occ[:,2]
    data['n20'] = occ[:,3]
    data['n21'] = occ[:,4]
    data['n22'] = occ[:,5]
    data['host_haloid'] = parents[:,1]
    print('Adding fields')
    rf.append_field(hlist,['lumtot', 'lum20', 'lbcg', 'ngals', 'n18',
                           'n19', 'n20', 'n21', 'n22', 'host_haloid'],
                    data, adtype)
    print('Writing file')
    fitsio.write('{0}/{1}_halos_{2}.fits'.format(outpath,prefix,bsize), hlist)

if __name__ == '__main__':
    
    basepath = '/nfs/slac/des/fs1/g/sims/jderose/l-addgals/catalogs/Buzzard/'
    halopaths = ['/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb1050/Octants_01/rockstar/halos_1050.fit',
                 '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb2600/Octants_01/rockstar/halos_2600.fit',
                 '/nfs/slac/g/ki/ki21/cosmo/mbusha/Simulations/Chinchillas-midway/Chinchilla-1/Lb4000/Octants_01/rockstar/halos_4000.fit']
    outpath = '/lustre/ki/pfs/jderose/l-addgals/catalogs/Buzzard/Catalog_l-addgals/'
    prefix = 'Buzzard'
    postfix = 'l-addgals'
    finalize_catalogs(basepath, prefix, postfix, outpath, halopaths)
    #occ = np.zeros((1808525,6))
    #lum = np.zeros((1808525,3))
    #combine_parents_list_fits(parentpath, halopath, prefix, outpath, bsize, occ, lum)

    
    

    
    
    
    
