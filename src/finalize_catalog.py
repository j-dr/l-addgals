from __future__ import print_function, division
import fitsio
from collections import deque
import numpy as np
import healpy as hp
from mpi4py import MPI
from scipy import spatial

def finalize_catalogs(basepath, prefix, postfix, outpath, halopaths, nh=9):

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    nio = size-nh-1
    srank = -1
    pixcat = None
    message = None

    tags = {'write':0, 'exit':1}

    if rank == 0:
        free = deque(range(nio,size))
        writing = deque()
        swaiting = deque()
        wwaiting = deque()
        done = []
        while True:
            #are we done?
            if len(done)==size: break

            #see if any write requests can be filled
            if len(wwaiting)>0:
                remove = []
                for w in wwaiting:
                    if w[1] not in writing:
                        comm.Send(tag=tags['write'], dest=w[0])
                        writing.append(w[1])
                        remove.append(w)
            
            for r in remove:
                wwaiting.remove(r)

            #Recieve requests
            status = MPI.Status()
            comm.Recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            if tag=tags['write']:
                wwaiting.append(status.Get_source())
            if tag=tags['exit']:
                done.append(status.Get_source())
                
            
    if (rank < nio) and (rank != 0):
        for i, bsize in enumerate(['1050', '2600', '4000']):
            #read in the halos for this box, sort by healpix cell number
            h = readHlist(hpath, usecols=[0,1,2,5,8,9,10], 
                          names=['id','pid','mvir','rvir','x','y','z'])
            occ = np.zeros(len(h))
            pixh = hp.vec2pix(2, h['x'], h['y'], h['z'])

            #index the halos by healpix
            pidx = pixh.argsort()
            pixh = pixh[pidx]
            h = h[pidx]
            pidx = pixh[1:]-pixh[:-1]
            pidx = np.where(pidx!=0)
            pidx = np.hstack([np.zeros(1),pidx])
            upix = pixh[pidx]
            pixmap = hp.ud_grade(np.arange(12*2**2),4)

            pixpaths = glob('{0}/Lb{1}_{2}/[0-9]*'.format(basepath, bsize, postfix))
            chunks = [pixpaths[i::nio] for i in range(nio)]

            for ppath in chunks[rank]:
                #build kdtree for this cell
                pix = int(ppath.split('/')[-1])
                hidx, = np.where(upix==pixmap[pix])
                hstart = pidx[hidx]
                if hstart==(len(pidx)-1):
                    hend = len(h)
                else:
                    hend = pidx[hidx+1]

                ht = spatial.KDTree(h[hstart:hend]['x','y','z'].view((h.dtype['x'],3)))

                pfiles = glob('{0}/*/hv_output/gal_ginfo1.dat'.format(ppath))
                
                for i, f in enumerate(pfiles):
                    pc, hdr = fitsio.read(f, header=True)
                    c = associate_halos(pc, h[hstart:hend], ht)
                    occ += c
                    
                    #request permission to write
                    status = MPI.Status()
                    comm.Send(pix, tag=tags['write'])
                    comm.Recv(tag=tags['write'], status=status)
                    write_fits(outpath, prefix, pix, pc, hdr)
                    
            comm.Send(tag=tags['exit'])


def associate_halos(galaxies, halos, tree):
    
    d, hid = tree.query(galaxies[['PX', 'PY', 'PZ']].view((galaxies.dtype['x'],3)))
    galaxies[galaxies['HALOID']!=-1]['RHALO'] = 0.0
    galaxies[galaxies['HALOID']==-1]['RHALO'] = 999999
    galaxies[galaxies['HALOID']==-1]['RHALO'] = d[galaxies['HALOID']==-1]
    galaxies[galaxies['HALOID']==-1]['HALOID'] = halos['id'][hid[galaxies['HALOID']==-1]]
    
    c, e = np.histogram(galaxies['HALOID'][galaxies['RHALO']<=halos[hid]], 
                        bins=np.arange(np.min(halos['id']),np.max(halos['id']))))

    return c

def write_fits(outpath, prefix, pix, catalog, header, nside=8):

    hpix = hp.vec2pix(nside, catalog['x'],catalog['y'],catalog['z'])
    upix = np.unique(hpix)

    for p in upix:
        fits = fitsio.FITS('{0}/{1}_{2}.fits'.format(outpath, prefix, p) ,'rw')
        try:
            pidx = hpix==p
            h = fits[-1].read_header()
            fits[-1].write_key('NAXIS2', h['NAXIS2']+len(catalog[pidx]))
            fits[-1].append(catalog[pidx])
        except:
            fits.write(catalog header=header)

    
    
    
    
