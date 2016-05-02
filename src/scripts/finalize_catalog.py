from __future__ import print_function, division
from SimulationAnalysis import readHlist
from collections import deque
from glob import glob
from mpi4py import MPI
import fast3tree as f3t
import numpy.lib.recfunctions as rf
import numpy as np
import pandas as pd
import healpy as hp
import fitsio
import warnings
import yaml
import time
import sys
import os


TZERO = None
def tprint(info, stime=None):
    global TZERO
    if TZERO is None:
        TZERO = time.time()

    tnow = time.time()
    print('[%8ds] %s' % (tnow-TZERO,info))
    
    return tnow



def finalize_catalogs(basepath, prefix, suffix, outpath, halopaths, ztrans, 
                      bzcut = [0.34, 0.9, 2.0], mmin=[3e12,3e12,2.4e13], zbuff=0.05, 
                      nside=8, justhalos=True, pixels=None, rmp_output=False,
                      lensing_output=False, skyfactory=True):

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    tags = {'write':0, 'fwrite':1, 'occ1050':2, 'occ2600':3, 'occ4000':4,
            'lum1050':5, 'lum2600':6, 'lum4000':7, 'exit':8, 'recv':9}
    bsizeenum = {'1050':0, '2600':1, '4000':2}
    message = None

    if rank == 0:
        writing = deque()
        wwaiting = deque()
        rwaiting = deque()
        done = []
        occ = [None, None, None]
        lum = [None, None, None]
        pixpaths = deque()
        for i, bsize in enumerate(['1050', '2600', '4000']):
            if skyfactory:
                pix = glob('{0}/Lb{1}_{2}/output/addgals/[0-9]*/0*'.format(basepath, bsize, suffix))
            else:
                pix = glob('{0}/Lb{1}_{2}/[0-9]*/0*'.format(basepath, bsize, suffix))
            zbins = np.unique(np.array([p.split('/')[-1] for p in pix]))
            if pixels!=None:
                pix = pixels
            else:
                pix = np.unique(np.array([p.split('/')[-2] for p in pix]))

            for zb in zbins:
                for p in pix:
                    pixpaths.append([bsize, p, zb])

        while True:
            #are we done?
            if len(done)==(size-1): break
            remove = []
            #see if any write requests can be filled
            if len(wwaiting)>0:
                for w in wwaiting:
                    tprint('    Rank {0} is waiting to write to pixels {1}'.format(w[0],w[1]))
                    write = True
                    for p in w[1]:
                        if p in writing: write=False
                    if write:
                        tprint('    Master gives rank {0} permission to write'.format(w[0]))
                        comm.send(message, tag=tags['write'], dest=w[0])
                        writing.extend(w[1])
                        remove.append(w)

            #see if any processors need another bin
            while len(rwaiting)>0:
                proc = rwaiting.popleft()
                if len(pixpaths)==0:
                    comm.send(None, tag=tags['recv'], dest=proc)
                    continue

                finfo = pixpaths.popleft()
                tprint('    Assigning rank {0} to bsize, pix, zbin: {1}, {2}, {3}'.format(proc, *finfo))
                comm.send(finfo, tag=tags['recv'], dest=proc)
            
            for r in remove:
                wwaiting.remove(r)

            #Need to process finished writing notifications more quickly 
            #so that requests don't build up
            while len(writing)>10:
                status = MPI.Status()
                tprint('    Master waiting for writing to finish...')
                message = comm.recv(source=MPI.ANY_SOURCE, tag=tags['fwrite'], status=status)
                for p in message:
                    writing.remove(p)

            #Recieve requests
            status = MPI.Status()
            tprint('    Master waiting to recieve messages')
            message = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            tprint('    Master recieved message')
            tag = status.Get_tag()
            if tag==tags['write']:
                wwaiting.append([status.Get_source(), message])
            elif tag==tags['fwrite']:
                for p in message:
                    writing.remove(p)
            elif tag==tags['recv']:
                rwaiting.append(status.Get_source())
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
            if i==0:
                zl = 0.0
            else:
                zl = bzcut[i-1]

            update_halo_file(halopaths[i], prefix, outpath, bsize, occ[i], 
                             lum[i], mmin[i], zl, bzcut[i], zbuff=zbuff)
            
    if rank != 0:
        bsize = None
        lastbsize = None
        occ = None
        lum = None
        while True:
            tprint('    {0}: Requesting new bin'.format(rank))
            comm.send(tag=tags['recv'])
            message = comm.recv(tag=tags['recv'])
            if message==None:
                tocc = tprint('    {0}: Sending occupations and finishing'.format(rank))
                if (occ!=None) and (lum!=None):
                    comm.send(occ, tag=tags['occ'+bsize])
                    comm.send(lum, tag=tags['lum'+bsize])
                    tprint('    {0}: Communication took {1}s'.format(rank, time.time()-tocc))
                break
            lastbsize = bsize
            bsize = message[0]
            pix = message[1]
            zbin = message[2]

            tstart = tprint('    {0}: Processing box, pix, zbin: {1}, {2}, {3}'.format(rank, bsize, pix, zbin))
            
            #If starting on a new box, read in the correct halo file
            if bsize!=lastbsize:
                if lastbsize!=None:
                    tocc = tprint('    {0}: Sending occupations'.format(rank))
                    comm.send(occ, tag=tags['occ'+lastbsize])
                    comm.send(lum, tag=tags['lum'+lastbsize])
                    tprint('    {0}: Communication took {1}s'.format(rank, time.time()-tocc))

                h = fitsio.read(halopaths[bsizeenum[bsize]], columns=['HALOID', 'MVIR', 'RVIR', 'HALOPX', 
                                                                      'HALOPY', 'HALOPZ', 'HOST_HALOID', 'Z'])

                if bsizeenum[bsize]==0:
                    hlz = 0.0
                else:
                    hlz = bzcut[bsizeenum[bsize]-1]

                #only want halos above mass cut in redshifts relevant for this box size
                h = h[h['MVIR']>mmin[bsizeenum[bsize]]]
                h = h[((hlz-zbuff)<=h['Z']) & (h['Z']<(bzcut[bsizeenum[bsize]]+zbuff))]

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
                tread = tprint('    {0}: Done reading and sorting halos'.format(rank))
                print('{0}: Reading took {1}s'.format(rank, tread-tstart))
            
            if skyfactory:
                f = '{0}/Lb{1}/output/addgals/{2}/{3}/{4}_{1}.{2}.{3}.fits'.format(basepath, bsize, pix, zbin, prefix)
            else:
                f = '{0}/Lb{1}_{2}/{3}/{4}/{5}_{1}.{3}.{4}.fits'.format(basepath, bsize, suffix, pix, zbin, prefix)

            #get redshift cutoffs for this bin
            ztidx, = np.where((ztrans['zbin']==zbin) & (ztrans['bsize']==bsize))
            assert(len(ztidx)==1)

            hidx, = np.where(upix==pixmap[int(pix)])
            hstart = pidx[hidx]

            if hidx==(len(pidx)-1):
                hend = len(h)
            else:
                hend = pidx[hidx+1]

            #build tree to associate galaxies and halos 
            with f3t.fast3tree(h[hstart:hend][['HALOPX','HALOPY','HALOPZ']].view((h.dtype['HALOPX'],3))) as ht:
                try:
                    pc, hdr = fitsio.read(f, header=True)
                except:
                    warnings.warn('Could not read file {0}, continuing...'.format(f))
                    continue

                #get galaxies in correct redshift bin
                pc = pc[((pc['PX']!=0) | (pc['PY']!=0) | (pc['PZ']!=0)) &
                        (ztrans['zmin'][ztidx]<=pc['Z']) & (pc['Z']<ztrans['zmax'][ztidx])]

                tsh = tprint('    {0}: Finding halos for pixel, redshift bin: {1}, {2}'.format(rank, pix, zbin))
                o, l = associate_halos(pc, h[hstart:hend], ht)
                tfh = tprint('    {0}: Finished finding halos for redshift bin. Took {1}s'.format(rank, time.time()-tsh))

                #update halo occupations
                occ[hstart:hend] += o
                lum[hstart:hend] += l
                
                #create unique ID
                pc['ID'] += int(bsize)*1000000000000 + int(pix)*10000000000 + int(zbin)*100000000

                #if writing catalogs, request permission to write
                if not justhalos:
                    hpix = hp.vec2pix(nside, pc['PX'],pc['PY'],pc['PZ'], nest=True)
                    up = np.unique(hpix)
                    status = MPI.Status()
                    treq = tprint('    {0}: Requesting permission to write to pix {1}'.format(rank, pix))
                    comm.send(up, tag=tags['write'])
                    message = comm.recv(tag=tags['write'], status=status)
                    
                    tprint('    {0}: Waited {1}s for permission'.format(rank, time.time()-treq))
                    twr = tprint('    {0}: Writing galaxies for pix, z bin: {1}, {2}'.format(rank, pix, zbin))
                    write_fits(outpath, prefix, pix, pc, hpix, up, hdr, 
                               rmp_output=rmp_output, lensing_output=lensing_output)
                    tdw = tprint('    {0}: Done writing galaxies for {1}, {2}. Took {3}s'.format(rank, pix, zbin, time.time()-twr))

                    comm.send(up, tag=tags['fwrite'])
                    tprint('    {0}: Message sent. Took {1}s'.format(rank, time.time()-tdw))

        comm.send(message, tag=tags['exit'])

def associate_halos(galaxies, halos, tree, rassoc=10):

    occ = np.zeros((len(halos),6))
    lum = np.zeros((len(halos),3))
    gpos = galaxies[['PX', 'PY', 'PZ']].view((galaxies.dtype['PX'],3))
    hpos = halos[['HALOPX','HALOPY','HALOPZ']].view((halos.dtype['HALOPX'],3))
    d = np.zeros(len(galaxies))
    hid = np.zeros(len(galaxies), dtype=np.int)
    
    #assign nearest halo to galaxy
    for i, gi in enumerate(gpos):
        try:
            gi = tree.query_radius(gpos[i,:],rassoc)
            di = np.sqrt((hpos[gi,0]-gpos[i,0])**2 + (hpos[gi,1]-gpos[i,1])**2 + (hpos[gi,2]-gpos[i,2])**2)
            d[i] = np.min(di)
            hid[i] = gi[di==d[i]]
        except:
            hid[i] = -1
            d[i] = -1
            
    #tprint('    Minimum distance to nearest halo: {0}'.format(np.min(d)))

    cen = galaxies['CENTRAL']==1
    nn = hid!=-1
    galaxies['RHALO'][cen] = 0.0
    galaxies['RHALO'][~cen] = 999999
    galaxies['RHALO'][~cen & nn] = d[~cen & nn]
    galaxies['HALOID'][~cen & nn] = halos['HALOID'][hid[~cen & nn]]
    galaxies['M200'][~cen & nn] = halos['MVIR'][hid[~cen & nn]]
    galaxies['M200'][cen & nn] = halos['MVIR'][hid[cen & nn]]
    bound = galaxies['RHALO']<=halos['RVIR'][hid]

    nmt, = np.where(halos['HALOID'][hid[cen]]==galaxies[cen]['HALOID'])
    tprint('    Number of centrals with unidentified halos: {0}'.format(len(galaxies[cen])-len(nmt)))

    galaxies[cen]['M200'] = halos[hid[cen]]['MVIR']
    lum[hid[cen],2] = galaxies[cen]['AMAG'][:,1]

    for i, mr in enumerate([-10,-18,-19,-20,-21,-22]):
        #get galaxies passing magnitude cut
        midx = galaxies['AMAG'][bound,1]<mr
        
        #sort halo indices of these galaxies
        hidx = hid[bound][midx].argsort()
        mhid = hid[bound][midx][hidx]
        if len(mhid)==0:continue

        #get magnitudes of these galaxies
        mgr = galaxies['AMAG'][bound,1][midx][hidx]

        #figure out which galaxies go with which halos
        hidx = mhid[1:]-mhid[:-1]
        hidx, = np.where(hidx!=0)
        hidx = np.hstack([np.zeros(1,dtype=np.int),hidx+1])
        uhid = mhid[hidx]
        hidx = np.hstack([hidx, np.array([len(mhid)])])
        if (mr==-10):
            for j, uid in enumerate(uhid):
                lum[uid,0]+=np.sum(mgr[hidx[j]:hidx[j+1]])
        elif (mr==-20):
            for j, uid in enumerate(uhid):
                lum[uid,1]+=np.sum(mgr[hidx[j]:hidx[j+1]])
        
        occ[uhid,i] += hidx[1:]-hidx[:-1]

    return occ, lum

def write_fits(outpath, prefix, pix, catalog, hpix, upix, header, nside=8, 
               rmp_output=False, lensing_output=False):

    for p in upix:
        fits = fitsio.FITS('{0}/truth/{1}.{2}.fits'.format(outpath, prefix, p) ,'rw')
        print('Writing to {0}/truth/{1}.{2}.fits'.format(outpath, prefix, p))
        pidx = hpix==p
        try:
            h = fits[-1].read_header()
            fits[-1].write_key('NAXIS2', h['NAXIS2']+len(catalog[pidx]))
            fits[-1].append(catalog[pidx])
        except:
            fits.write(catalog[pidx], header=header)

        fits.close()

        if rmp_output:
            columns = ['ID', 'RA', 'DEC', 'OMAG', 'OMAGERR', 'M200', 'CENTRAL', 'RHALO', 'Z']
            fits = fitsio.FITS('{0}/rmp/{1}.{2}.rmp.fits'.format(outpath, prefix, p) ,'rw')
            try:
                h = fits[-1].read_header()
                fits[-1].write_key('NAXIS2', h['NAXIS2']+len(catalog[pidx]))
                fits[-1].append(catalog[columns][pidx])
            except:
                fits.write(catalog[columns][pidx], header=header)
                
            fits.close()

        if lensing_output:
            columns = ['ID', 'PX', 'PY', 'PZ']
            fits = fitsio.FITS('{0}/lens/{1}.{2}.lens.fits'.format(outpath, prefix, p) ,'rw')
            try:
                h = fits[-1].read_header()
                fits[-1].write_key('NAXIS2', h['NAXIS2']+len(catalog[pidx]))
                fits[-1].append(catalog[columns][pidx])
            except:
                fits.write(catalog[columns][pidx], header=header)
                
            fits.close()


def update_halo_file(halopath, prefix, outpath, bsize, occ, lum, mmin, zmin, zmax, zbuff=0.05):

    h, header = fitsio.read(halopath, header=True)

    h = h[h['MVIR']>mmin]
    h = h[((zmin-zbuff)<=h['Z']) & (h['Z']<(zmax+zbuff))]
    pixh = hp.vec2pix(2, h['HALOPX'], h['HALOPY'], h['HALOPZ'])
    pidx = pixh.argsort()
    pixh = pixh[pidx]
    h = h[pidx]
    h['LUMTOT'] = lum[:,0]
    h['LUM20'] = lum[:,1]
    h['LBCG'] = lum[:,2]
    h['NGALS'] = occ[:,0]
    h['N18'] = occ[:,1]
    h['N19'] = occ[:,2]
    h['N20'] = occ[:,3]
    h['N21'] = occ[:,4]
    h['N22'] = occ[:,5]

    pidx = pixh[1:]-pixh[:-1]
    pidx, = np.where(pidx!=0)
    pidx = np.hstack([np.zeros(1,dtype=np.int),pidx+1])
    upix = pixh[pidx]
    
    for i, p in enumerate(upix):
        start = pidx[i]
        if i==len(upix)-1:
            end = pidx[i+1]
        else:
            end = len(h)

        fits = fitsio.FITS('{0}/{1}_halos.{2}.fits'.format(outpath, prefix, p) ,'rw')
        print('Writing to {0}/{1}_halos.{2}.fits'.format(outpath, prefix, p))
        try:
            hdr = fits[-1].read_header()
            fits[-1].write_key('NAXIS2', hdr['NAXIS2']+len(h[start:end]))
            fits[-1].append(h[start:end])
        except:
            fits.write(h[start:end], header=header)

        fits.close()

        if lensing_output:
            columns = ['HALOID', 'HALOPX', 'HALOPY', 'HALOPZ']
            fits = fitsio.FITS('{0}/lens/{1}.{2}_halos.lens.fits'.format(outpath, prefix, p) ,'rw')
            hp = h[columns][start:end]
            hp.dtype.names = 'ID', 'PX', 'PY', 'PZ'
            try:
                hdr = fits[-1].read_header()
                fits[-1].write_key('NAXIS2', hdr['NAXIS2']+len(h[start:end]))
                fits[-1].append(hp)
            except:
                fits.write(hp)
                
            fits.close()



def join_halofiles(basepath, mmin=5e12):

    lnames = ['ID', 'MVIR', 'VMAX', 'VRMS', 'VRMS', 'RVIR', 'RS', 'NP', 'HALOPX', 'HALOPY', 'HALOPZ',
              'HALOVX', 'HALOVY', 'HALOVZ', 'JX', 'JY', 'JZ', 'SPIN', 'RS_KLYPIN', 'MVIR_ALL', 
              'M200b', 'M200c', 'M500c', 'M2500c', 'XOFF', 'VOFF', 'LAMBDA', 'B_TO_A', 'C_TO_A',
              'AX', 'AY', 'AZ', 'VIRIAL_RATIO']

    pnames = ['ID','MVIR','PID']
    
    lusecols = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
    pusecols = [0,2,41]

    print('Reading parents')
    parents = pd.read_csv("{0}/cut_reform_out_0.parents".format(basepath), usecols=pusecols, names = lnames, comment='#', sep=' ')
    parents = parents.to_records(index=False)
    parents = parents[parents['MVIR']>mmin]

    print('Reading list')
    hlist = pd.read_csv("{0}/cut_reform_out_0.list".format(basepath), usecols=lusecols, names=pnames, comment='#', sep=' ')
    hlist = hlist.to_records(index=False)
    hlist = hlist[hlist['MVIR']>mmin]

    print('Joining files')
    hlist = rf.join_by('ID', hlist, parents[['ID', 'PID']], r1postfix=None, r2postfix=None, usemask=False, asrecarray=True)
    adtype = np.dtype([('LUMTOT',np.float), ('LUM20',np.float), ('LBCG', np.float),
                       ('NGALS',np.int), ('N18',np.int), ('N19',np.int), ('N20',np.int),
                       ('N21',np.int), ('N22',np.int)])

    print('Adding fields')
    rf.append_field(hlist,['LUMTOT', 'LUM20', 'LBCG', 'NGALS', 'N18',
                           'N19', 'N20', 'N21', 'N22'], dtype=adtype)

    print('Writing file')
    fitsio.write('{0}/out_0.fits'.format(basepath), hlist)


if __name__ == '__main__':
    
    cfgfile = sys.argv[1]
    with open(cfgfile, 'r') as fp:
        cfg = yaml.load(fp)

    outpath = cfg['outpath']
    halopaths = cfg['halopaths']
    basepath = cfg['basepath']
    prefix = cfg['prefix']
    suffix = cfg['suffix']
    zbinfile = cfg['zbinfile']
    zbins = fitsio.read(zbinfile)
    if 'pixels' in cfg.keys():
        pixels = cfg['pixels']
    else:
        pixels = None

    if 'justhalos' in cfg.keys():
        justhalos = cfg['justhalos']
    else:
        justhalos = False

    if 'rmp_output' in cfg.keys():
        rmp_output = True
    else:
        rmp_output = False

    if 'lensing_output' in cfg.keys():
        lensing_output = True
    else:
        lensing_output = False
        
    if 'skyfactory' in cfg.keys():
        skyfactory = cfg['skyfactory']
    else:
        skyfactory = False

    if 'mmin' in cfg.keys():
        mmin = [float(m) for m in mmin]
    else:
        mmin = [3e12,3e12,2.4e13]

    try:
        os.makedirs(outpath)
    except:
        pass

    jhalopaths = []
    for i, hp in enumerate(halopaths):
        hs = hp.split('.')
        if 'fit' not in hs[-1]:
            jhalopaths.append(join_halofiles('.'.join(hs[:-1]), mmin=mmin[i]))
        else:
            jhalopaths.append(hp)
            
    finalize_catalogs(basepath, prefix, suffix, outpath, jhalopaths, zbins, 
                      justhalos=justhalos, pixels=pixels, lensing_output=lensing_output,
                      rmp_output=rmp_output, skyfactory=skyfactory)
