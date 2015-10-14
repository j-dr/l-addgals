import numpy as np
import pyfits
import os
import sys

cspeed = 2.99792e5

nbins = 10
rmin = 0.1
rmax = 50

def count_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i+1

# short function for opening a fits file and getting the data
def read_fits(fname):
    hdu = pyfits.open(fname)
    d = hdu[1].data

    return d

def run_wp(nd1, nd2, nr1, nr2, ddfile, drfile, rrfile, wpfile):

    exe = '/afs/slac.stanford.edu/u/ki/mbusha/projects/berlind_wp/wprp' 

    cmd = exe+' '+str(nbins)+' '+str(nd1)+' '+str(nd2)+' '+str(nr1)+' '+str(nr2)+' '+ddfile+' '+drfile+' '+drfile+' '+rrfile+' > '+wpfile
    print cmd
    sys.stdout.flush()
    os.system(cmd)

def run_ddrpi(infile1, infile2, outfile):

    exe1 = '/afs/slac.stanford.edu/u/ki/mbusha/projects/berlind_wp/DDrppi'

    cmd1 = exe1+' '+str(rmin)+' '+str(rmax)+' '+str(nbins)+' '+infile1+' '+infile2+' > '+outfile
    print cmd1
    sys.stdout.flush()
    os.system(cmd1)

# class defining a galaxy object
class Galaxy(object):

    def __init__(self):

          #self.name = name

          return

    #read the data
    def read_data(self, halofile):

        print halofile
        d = read_fits(halofile)

        self.ra = d['ra'][:]
        self.dec = d['dec'][:]
        self.z = d['z'][:]

        self.mr = d['amag'][:,2]
        self.gmr = d['amag'][:,1] - d['amag'][:,2]

        return

    # determines if the object is red or blue
    def is_red(self, i):

        if (0.21 - self.gmr[i] < 0.03*self.mr[i]):

            return 1

        return 0

    def write_wp_info(self, maglim, zlim, outall, outred, outblue):

        # save all, red, and blue objects in separate files
        fall = open(outall, 'w')
        fred = open(outred, 'w')
        fblue = open(outblue, 'w')

        nObj = len(self.ra)
        print 'nObj: {0}'.format(nObj)
        for i in range(nObj):

            if (self.z[i] > zlim):
                continue
            if (self.mr[i] > maglim):
                continue

            s = str(self.ra[i]) + ' '+str(self.dec[i]) + ' ' + str(cspeed*self.z[i])+'\n'
            fall.write(s)
            if self.is_red(i):
                fred.write(s)
            else:
                fblue.write(s)

        fall.close()
        fblue.close()
        fred.close()

        return 

####################################################
# start our main program
####################################################

# fbase = str(sys.argv[1])
# inpath = str(sys.argv[2])
# outpath = str(sys.argv[3])
Mr = int(sys.argv[1])

fbase = 'galaxies_wp_cuts.fit'
inpath= '/nfs/slac/g/ki/ki23/des/jderose/addgals/catalogs/chinchilla-tuning/SVA/validation/wp_cat/'
outpath= '/nfs/slac/g/ki/ki23/des/jderose/addgals/catalogs/chinchilla-tuning/SVA/validation/wp_cat/'
#inpath = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla-1/Catalog_v2.0/wp_cat/'
#outpath = '/nfs/slac/g/ki/ki23/des/mbusha/catalogs/Chinchilla-1/Catalog_v2.0/wp_cat/'
#Mr = 21

if not os.path.exists(outpath):
    os.makedirs(outpath)
    

# define our output files
tmppath = '/nfs/slac/g/ki/ki22/cosmo/jderose/scratch/'
infile = inpath+fbase
refbase = tmppath+'/reformatted_'+fbase+'.dat'
ddbase = tmppath+'/DDrpi_'+fbase+'.dat'
drbase = tmppath+'/DRrpi_'+fbase+'.dat'
wpbase = outpath+'/wp_'+fbase

# location of our random files
ranbase = '/nfs/slac/g/ki/ki23/des/jderose/addgals/catalogs/chinchilla-tuning/SVA/validation/random/ranfile_'
rrbase = '/nfs/slac/g/ki/ki23/des/jderose/addgals/catalogs/chinchilla-tuning/SVA/validation/random/rrpi'

# definitions for our volume limited samples
# Mr = [18, 19, 20, 21]
# zmax = [12500.0/cspeed, 19250.0/cspeed, 31900.0/cspeed, 47650./cspeed]
nr = 1000000
if Mr == 18:
    zmax = 12500.0/cspeed
elif Mr == 19:
    zmax = 19250.0/cspeed
elif Mr == 20:
    zmax = 31900.0/cspeed
elif Mr == 21:
    zmax = 47650.0/cspeed
else:
    print "Error!  You must specify Mr = 18, 19, 20, or 21.  You supplied", Mr
    sys.exit(0)

print "Creating structure..."
g = Galaxy()

print "Reading the data..."
print infile
g.read_data(infile)

print "Doing the calculation..."

# define the files we're saving the reformatted information to
ref_all = refbase+str(Mr)+'_all.dat'
ref_red = refbase+str(Mr)+'_red.dat'
ref_blue = refbase+str(Mr)+'_blue.dat'

# reformat the files
g.write_wp_info(-1*Mr, zmax, ref_all, ref_red, ref_blue)
nd_all = count_lines(ref_all)
#nd_red = count_lines(ref_red)
#nd_blue = count_lines(ref_blue)

# where we'll save the DD and DR files
dd_all = ddbase+str(Mr)+'_all.dat'
dd_red = ddbase+str(Mr)+'_red.dat'
dd_blue = ddbase+str(Mr)+'_blue.dat'

dr_all = drbase+str(Mr)+'_all.dat'
dr_red = drbase+str(Mr)+'_red.dat'
dr_blue = drbase+str(Mr)+'_blue.dat'

# crate the dd files
print "Creating the DDrpi files..."
run_ddrpi(ref_all, ref_all, dd_all)
#run_ddrpi(ref_red, ref_red, dd_red)
#run_ddrpi(ref_blue, ref_blue, dd_blue)

# the paths to the random files
ranfile = ranbase+str(Mr)+'.dat'
rrfile = rrbase+str(Mr)+'.dat'

# create the dr files
print "Creating the DRrpi files..."
run_ddrpi(ref_all, ranfile, dr_all)
#run_ddrpi(ref_red, ranfile, dr_red)
#run_ddrpi(ref_blue, ranfile, dr_blue)

# create the wp command
print "Creating the wp files..."
wp_all = wpbase+'_'+str(Mr)+'_all.dat'
wp_red = wpbase+'_'+str(Mr)+'_red.dat'
wp_blue = wpbase+'_'+str(Mr)+'_blue.dat'
run_wp(nd_all, nd_all, nr, nr, dd_all, dr_all, rrfile, wp_all)
#run_wp(nd_red, nd_red, nr, nr, dd_red, dr_red, rrfile, wp_red)
#run_wp(nd_blue, nd_blue, nr, nr, dd_blue, dr_blue, rrfile, wp_blue)

