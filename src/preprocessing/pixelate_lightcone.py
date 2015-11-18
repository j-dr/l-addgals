from __future__ import print_function, division
from collections import namedtuple
import numpy as np
import healpy as hp
import struct

__GadgetHeader_fmt = '6I6dddii6Iiiddddii6Ii'
__finenside = 8192

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
            print( header )
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

def write_to_redshift_hpix_cells(filepath, outbase, filenside=64, rmin=0, rmax=4000, rstep=25):
    """
    Read in gadget particle block, and write to the correct healpix/redshift
    cell files
    """

    rbins = np.linspace(rmin,rmax,(rmax-rmin)//rstep)
    
    hdr, pos, vel, ids = readGadgetSnapshot(filepath, read_pos=True, read_vel=True,\
                                                read_id=True, lgadget=True)

    r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
    pix = hp.vec2pix(filenside, pos[:,0], pos[:,1], pos[:,2], nest=True)
    upix = np.sort(np.unique(pix))
    
    idxdtype = np.dtype([('pidx', np.int32), ('ridx', np.int32)])
    bins = np.ndarray(len(pos), dtype=idxdtype)

    #Index by pixel and radial bin
    bins['ridx'] = np.digitize(r, rbins)
    bins['pidx'] = np.digitize(pix, upix)

    #sort indices
    idx = bins.argsort(order=['pidx', 'ridx'])

    pos = pos[idx,:]
    vel = vel[idx,:]
    ids = ids[idx]
    pix = pix[idx]
    bins = bins[idx]
    del idx

    #create index into pixel/z bins
    rinc = bins['ridx'][1:]-bins['ridx'][:-1]
    pinc = bins['pidx'][1:]-bins['pidx'][:-1]

    idx, = np.where((rinc!=0) | (pinc!=0))
    npart = idx[1:]-idx[:-1]

    #Write to disk
    for i, start in enumerate(idx):
        if i==(len(idx)-1):
            end = len(bins)-1
        else:
            end = idx[i+1]

        deltan = end-start

        with open(outbase+'_{0}_{1}.pos'.format(*bins[start]), 'w+b') as fp:
            #keep track of number of particles in bin
            try:
                npt = np.fromstring(fp.read(8), np.int64)[0]+deltan
            except:
                npt = deltan
            fp.seek(0)
            fp.write(struct.pack('q',npt))

            #append particles to end of file
            fp.seek(0,2)
            fmt = 'f'*3*deltan
            bin = struct.pack(fmt, *pos[start:end,:].flatten())
            fp.write(bin)
            
        with open(outbase+'_{0}_{1}.vel'.format(*bins[start]), 'ab') as fp:
            fmt = 'f'*3*deltan
            bin = struct.pack(fmt, *vel[start:end,:].flatten())
            fp.write(bin)
        
        with open(outbase+'_{0}_{1}.id'.format(*bins[start]), 'ab') as fp:
            fmt = 'q'*deltan
            bin = struct.pack(fmt, *ids[start:end])
            fp.write(bin)


def read_unprocessed_cell(filebase, pix, zbin):
    
    with open(filebase+'_{0}_{1}.pos'.format(pix,zbin), 'rb') as fp:
        fmt = np.dtype(np.float32)
        npart = np.fromstring(fp.read(8), np.int64)[0]
        pos = np.fromstring(fp.read(npart*3*fmt.itemsize), fmt).reshape((npart,3))

    with open(filebase+'_{0}_{1}.vel'.format(pix,zbin), 'rb') as fp:
        fmt = np.dtype(np.float32)
        vel = np.fromstring(fp.read(npart*3*fmt.itemsize), fmt).reshape((npart,3))

    with open(filebase+'_{0}_{1}.id'.format(pix,zbin), 'rb') as fp:
        fmt = np.dtype(np.int64)
        ids = np.fromstring(fp.read(npart*fmt.itemsize), fmt)

    return npart, pos, vel, ids


def read_processed_cell(filebase, zbin, filenside=64, read_pos=True, \
                            read_vel=True, read_ids=True):

    hdrfmt = 'fqf'
    idxfmt = np.dtype(np.int64)
    to_read = [read_pos, read_vel, read_ids]
    fmt = [np.dtype(np.float32), np.dtype(np.float32), np.dtype(np.int64)]
    item_per_row = [3,3,1]
    data = []
    filenpix  = 12*filenside**2

    with open(filebase+'_{0}'.format(zbin), 'rb') as fp:
        #read the header
        h = list(struct.unpack(hdrfmt, \
                fp.read(struct.calcsize(hdrfmt))))
        npart = h[1]
        data.append(h)
        #read the peano index
        idx = np.fromstring(fp.read(8*filenpix), idxfmt)
        data.append(idx)

        for i, r in enumerate(to_read):
            if r:
                data.append(np.fromstring(fp.read(npart*item_per_row[i]*fmt[i].itemsize), fmt[i]))

    return data


def nest2peano(pix, order):
  subpix = np.array([ [ 0, 1, 3, 2 ], [ 3, 0, 2, 1 ], [ 2, 3, 1, 0 ], [ 1, 2, 0, 3 ],\
                      [ 0, 3, 1, 2 ], [ 1, 0, 2, 3 ], [ 2, 1, 3, 0 ], [ 3, 2, 0, 1 ] ])
  subpath = np.array([ [ 4, 0, 6, 0 ], [ 7, 5, 1, 1 ], [ 2, 4, 2, 6 ], [ 3, 3, 7, 5 ],\
                       [ 0, 2, 4, 4 ], [ 5, 1, 5, 3 ], [ 6, 6, 0, 2 ], [ 1, 7, 3, 7 ] ])
  face2path = np.array([ 2, 5, 2, 5, 3, 6, 3, 6, 2, 3, 2, 3 ])
  face2peanoface = np.array([ 0, 5, 6, 11, 10, 1, 4, 7, 2, 3, 8, 9 ])
  
  npix_ = 12*(1 << (2*order))
  assert((pix >= 0).all() and (pix < npix_).all())
  
  face = pix>>(2*order)
  path = face2path[face]
  result = np.zeros(len(pix), dtype=np.int64);
  shifts = np.arange(0, 2*order-1, 2)

  for shift in shifts[::-1]:
      spix = (pix>>shift) & 0x3
      result <<= 2
      result |= subpix[path,spix]
      path = subpath[path,spix]

  return result + ((face2peanoface[face])<<(2*order));


def process_redshift_cell(filebase, zbin, header, filenside=64):

    partnside = 4096
    partorder = int(np.log2(partnside))
    fileorder = int(np.log2(filenside))
    coarsenpix = 12*filenside**2
    coarsepix = np.arange(coarsenpix)
    
    #order pixels by peano indices
    peano = nest2peano(coarsepix, fileorder)
    pidx = peano.argsort()
    coarsepix = coarsepix[pidx]

    peano = peano[pidx]
    idx = np.zeros(12*filenside**2, dtype=np.int64)

    hdrfmt = 'fqf'
    npart = 0
    with open(filebase+'_{0}'.format(zbin), 'wb') as fp:
        fp.write(struct.pack(hdrfmt, header['L_b'], npart, header['M_p']))
        fp.write(struct.pack('q'*coarsenpix, *idx))

        for i, hpix in enumerate(coarsepix):
            try:
                npart, pos, vel, ids = read_unprocessed_cell(filebase, hpix, zbin)
            except:
                continue

            idx[i] = npart
            pix = hp.vec2pix(partnside, pos[:,0], pos[:,1], pos[:,2], nest=True)
            peano = nest2peano(pix, partorder)
            del pix
            
            pidx = peano.argsort()
            pos = pos[pidx,:]
            vel = vel[pidx,:]
            ids = ids[pidx]
            del pidx, peano
            
            #write sorted particles
            fp.write(struct.pack('f'*3*npart, *pos.flatten()))
            #write sorted velocities
            fp.write(struct.pack('f'*3*npart, *vel.flatten()))
            #write sorted inds
            fp.write(struct.pack('q'*npart, *ids))

        fp.seek(0)
        fp.write(struct.pack(hdrfmt, header['L_b'], np.sum(idx), header['M_p']))
        fp.write(struct.pack('q'*coarsenpix, *idx))

