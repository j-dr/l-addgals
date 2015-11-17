from __future__ import print_function, division
from collections import namedtuple
import numpy as np
import healpy as hp
import struct

__GadgetHeader_fmt = '6I6dddii6Iiiddddii6Ii'

GadgetHeader = namedtuple('GadgetHeader', \
        'npart mass time redshift flag_sfr flag_feedback npartTotal flag_cooling num_files BoxSize Omega0 OmegaLambda HubbleParam flag_age flag_metals NallHW flag_entr_ics')

def readGadgetSnapshot(filename, read_pos=False, read_vel=False, read_id=False,\
        read_mass=False, print_header=False, single_type=-1, lgadget=False):
    """
    This function reads the Gadget-2 snapshot file.

    Parameters
    ----------
    filename : str
        path to the input file
    read_pos : bool, optional
        Whether to read the positions or not. Default is false.
    read_vel : bool, optional
        Whether to read the velocities or not. Default is false.
    read_id : bool, optional
        Whether to read the particle IDs or not. Default is false.
    read_mass : bool, optional
        Whether to read the masses or not. Default is false.
    print_header : bool, optional
        Whether to print out the header or not. Default is false.
    single_type : int, optional
        Set to -1 (default) to read in all particle types. 
        Set to 0--5 to read in only the corresponding particle type.
    lgadget : bool, optional
        Set to True if the particle file comes from l-gadget. 
        Default is false.

    Returns
    -------
    ret : tuple
        A tuple of the requested data. 
        The first item in the returned tuple is always the header.
        The header is in the GadgetHeader namedtuple format.
    """
    blocks_to_read = (read_pos, read_vel, read_id, read_mass)
    ret = []
    with open(filename, 'rb') as f:
        f.seek(4, 1)
        h = list(struct.unpack(__GadgetHeader_fmt, \
                f.read(struct.calcsize(__GadgetHeader_fmt))))
        if lgadget:
            h[30] = 0
            h[31] = h[18]
            h[18] = 0
            single_type = 1
        h = tuple(h)
        header = GadgetHeader._make((h[0:6],) + (h[6:12],) + h[12:16] \
                + (h[16:22],) + h[22:30] + (h[30:36],) + h[36:])
        if print_header:
            print header
        if not any(blocks_to_read):
            return header
        ret.append(header)
        f.seek(256 - struct.calcsize(__GadgetHeader_fmt), 1)
        f.seek(4, 1)
        #
        mass_npart = [0 if m else n for m, n in zip(header.mass, header.npart)]
        if single_type not in range(6):
            single_type = -1
        #
        for i, b in enumerate(blocks_to_read):
            if i < 2:
                fmt = np.dtype(np.float32)
                item_per_part = 3
                npart = header.npart
            elif i==2:
                fmt = np.dtype(np.uint64) if lgadget or any(header.NallHW) \
                        else np.dtype(np.uint32)
                item_per_part = 1
                npart = header.npart
            elif i==3:
                fmt = np.dtype(np.float32)
                if sum(mass_npart) == 0:
                    ret.append(np.array([], fmt))
                    break
                item_per_part = 1
                npart = mass_npart
            size_per_part = item_per_part*fmt.itemsize
            #
            f.seek(4, 1)
            if not b:
                f.seek(sum(npart)*size_per_part, 1)
            else:
                if single_type > -1:
                    f.seek(sum(npart[:single_type])*size_per_part, 1)
                    npart_this = npart[single_type]
                else:
                    npart_this = sum(npart)
                data = np.fromstring(f.read(npart_this*size_per_part), fmt)
                if item_per_part > 1:
                    data.shape = (npart_this, item_per_part)
                ret.append(data)
                if not any(blocks_to_read[i+1:]):
                    break
                if single_type > -1:
                    f.seek(sum(npart[single_type+1:])*size_per_part, 1)
            f.seek(4, 1)
    #
    return tuple(ret)

if __name__=='__main__':

    filenside = 64
    outbase = 'pixelated_lightcone/lightcone_snapshot_'
    rbins = np.linspace(0,4000,4000//25)
    
    hdr, pos, vel, ids = readGadgetSnapshot(filepath, read_pos=True, read_vel=True, read_id=True, lgadget=True)

    r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
    pix = hp.vec2pix(filenside, pos[:,0], pos[:,1], pos[:,2])
    upix = np.sort(np.unique(pix))
    
    idxdtype = np.dtype([('pidx', np.int32), ('ridx', np.int32)])
    bins = np.ndarray(len(pos), dtype=idxdtype)

    bins['ridx'] = np.digitize(r, rbins)
    bins['pidx'] = np.digitize(pix, upix)

    idx = bins.argsort(order=['pidx', 'ridx'])

    pos = pos[idx,:]
    vel = vel[idx,:]
    ids = ids[idx]
    pix = pix[idx]
    bins = bins[idx]
    del idx

    pb = bins[0]['pidx']    
    rb = bins[0]['ridx']
    idx = [0]
    fbins = [(pb,rb)]

    for i, b in enumerate(bins):
        if (b[0]!=pb) | (b[1]!=rb):
            rb = b['ridx']
            pb = b['pidx']
            idx.append(i)
            fbins.append((pb,rb))

    idx = np.array(idx)

    for i, start in enumerate(idx):
        if i==(len(idx)-1):
            end = len(bins)-1
        else:
            end = idx[i+1]

        with open(lightconebase+'_{0}_{1}.pos'.format(*fbins[i]), 'wb') as fp:
            fmt = 'f'*3*(end-start)
            bin = struct.pack(fmt, *pos[start:end,:].flatten())
            fp.write(bin)
            
        with open(lightconebase+'_{0}_{1}.vel'.format(*fbins[i]), 'wb') as fp:
            fmt = 'f'*3*(end-start)
            bin = struct.pack(fmt, *vel[start:end,:].flatten())
            fp.write(bin)
        
        with open(lightconebase+'_{0}_{1}.id'.format(*fbins[i]), 'wb') as fp:
            fmt = 'q'(end-start)
            bin = struct.pack(fmt, *ids[start:end])
            fp.write(bin)

        
            
        
    
    
