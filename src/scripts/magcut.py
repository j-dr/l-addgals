#!/usr/bin/env python
from __future__ import print_function, division
import fitsio
import csv
import numpy as np
from glob import glob
import argparse

def cutmag(path, cmock, band, minmag):
    
    #Get all files to cut
    hv_path = '/'.join(path.split('/')[:-1])+'/hv_output'
    ff = glob(path+'/*.fit')
    #ff = [path+'/'+f for f in ff]
    df = glob(hv_path+'/*info*.dat')
    #df = [hv_path+'/'+f for f in df]

    #Select mock catalog whose mags we will use to make the cut
    cf = [f for f in ff if cmock in f][0]
    m = fitsio.read(cf, ext=1)
    
    ii = np.where((m['OMAG'][:,band]<minmag) | (m['OMAG'][:,band]!=m['OMAG'][:,band]))
    mc = m[ii]
    fitsio.write(cf+'.cut', mc)

    #Cut all other files 
    for f in ff:
        if f==cf: continue
        m = fitsio.read(f, ext=1)
        mc = m[ii]
        fitsio.write(f+'.cut', mc)

    for f in df:
        m = np.genfromtxt(f)
        mc = m[ii]
        with open(f+'.cut', 'w') as fp:
            w = csv.writer(fp, delimiter='\t')
            w.writerows(mc)


if __name__=='__main__':
    p = argparse.ArgumentParser()
    p.add_argument('path', type=str, help='path to .fit files')
    p.add_argument('cmock', type=str, help='name of mock to use for cut')
    p.add_argument('band', type=int, help='index of band to use for cut')
    p.add_argument('minmag', type=float, help='magnitude cut')

    args = p.parse_args()
    cutmag(args.path, args.cmock, args.band, args.minmag)
        
    
    

    

    

    
