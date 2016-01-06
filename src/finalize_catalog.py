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

    tags = {'send':1, 'recv':2, 'write':3, 'exit':4}

    if rank == 0:
        free = deque(range(nio,size))
        writing = deque()
        swaiting = deque()
        wwaiting = deque()
        done = []
        while True:
            #are we done?
            if len(done)==size:
                for r in range(nio,size):
                    comm.send('exit',tag=tags['recv'])
                break

            #see if any halo requests can be filled
            if (len(swaiting)>0) and (len(free)>0):
                comm.Send(free.popleft(), tag=tags['send'], dest=swaiting.popleft())

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
            if tag=tags['send']:
                swaiting.append(status.Get_source())
            if tag=tags['write']:
                wwaiting.append(status.Get_source())
            if tag=tags['exit']:
                done.append(status.Get_source())
                
            
    if (rank < nio) and (rank != 0):
        for i, bsize in enumerate(['1050', '2600', '4000']):
            pixpaths = glob('{0}/Lb{1}_{2}/[0-9]*'.format(basepath, bsize, postfix))
            
            chunks = [cfiles[i::nio] for i in range(nio)]
            
            for ppath in chunks[rank]:
                pix = ppath.split('/')[-1]
                pfiles = glob('{0}/*/hv_output/gal_ginfo1.dat'.format(ppath))
                
                for i, f in enumerate(pfiles):
                    if i==0:
                        pixcat, hdr = fitsio.read(f, header=True)
                    else:
                        pc = fitsio.read(f)
                        pixcat = np.hstack([pixcat, pc])

                status = MPI.Status()
                comm.Send(tag=tags['send']) #tell queen catalog ready for halo association
                comm.Recv(srank, source=0, tag=tags['send'], status=status) #wait for permission to send
                comm.send([hdr,pixcat,bsize], dest=srank, tag=tags['recv'])

    if rank >= nio:
        halos = []
        trees = []
        for hpath in halopaths:
            h = readHlist(hpath, usecols=[0,1,2,5,8,9,10], 
                          names=['id','pid','mvir','rvir','x','y','z'])
            ht = spatial.KDTree(h['x','y','z'].view((h.dtype['x'],3)))
            halos.append(h)
            trees.append(ht)

        while True:
            comm.recv(pixcat, source=MPI.ANY_SOURCE, tag=tags['recv'], status=status)
            if pixcat=='exit':
                break

            #find halos
            associate_halos(pixcat[1], halos[pixcat[2]], trees[pixcat[2]])
            #request permission to write
            status = MPI.Status()
            comm.Send(pix, tag=tags['write'])
            comm.Recv(tag=tags['write'], status=status)
            
            fits = fitsio.FITS('{0}/{1}_{2}.fits'.format(outpath, prefix, pix) ,'rw')
            try:
                fits[-1].append(pixcat[1])
            except:
                fits.write(pixcat[1], header=pixcat[0])

def associate_halos(galaxies, halos, tree):
    
    d, hid = tree.query(galaxies[['PX', 'PY', 'PZ']].view((galaxies.dtype['x'],3)))
    galaxies[galaxies['HALOID']==-1]['RHALO'] = d[galaxies['HALOID']==-1]
    galaxies[galaxies['HALOID']==-1]['HALOID'] = hid[galaxies['HALOID']==-1]

    
    
    
