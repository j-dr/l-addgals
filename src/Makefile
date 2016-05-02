#compiler
CXX=g++
CC=${CXX}
OPT = -O3 -Wno-deprecated -ffast-math -g #-ffast-math
WARN =
#END_FLAG = -DBIGENDIAN
END_FLAG = -DLITTLEENDIAN


NR_LIB = -L/home/jderose/src/recipes
NR_INC = -I/home/jderose/src/recipes
KCORRECT_LIB = -L${KCORRECT_DIR}/lib
KCORRECT_INC = -I${KCORRECT_DIR}/include
CCFITS_LIB   = -L/home/jderose/usr/local/lib
CCFITS_INC   = -I/home/jderose/usr/local/include

RM = rm -f
LD = /usr/bin/ld

#directories to look for include files
INCL = -I. -I${LOCAL_INC} ${NR_INC} ${KCORRECT_INC} ${CCFITS_INC}

#directories to look for libraries
LDFLAGS = -L${LOCAL_LIB} ${NR_LIB} ${KCORRECT_LIB} ${CCFITS_LIB}

#libraries to link
LOADLIBES = -lANN -lrecipes -lcfitsio -lkcorrect -lCCfits -lgsl -lgslcblas

#-------------> for Catalog Production <------------#
#catalog_type="classic"
catalog_type="BCC"
#catalog_type="SHAM_TEST"
#catalog_type="snapshot"
#catalog_type="just_colors"
#catalog_type="snapshot_colors"

ifeq ($(catalog_type),"classic")
CPPFLAGS = $(WARN) $(OPT) $(INCL) -DCATALOG $(END_FLAG)
CPPFLAGS += -DCOLORS            #Do we add colors to the galaxies?
CPPFLAGS += -DBCGS              #Add BCGs in accordance with Sarah's observations
CPPFLAGS += -DCOLORS_FROM_RELATIVE_DENSITY #Match SEDs using relative galaxy density in zbins instaed of absoloute densities
#CPPFLAGS += -DDENSPDF_IN_MAG_BIN  #Use DENSPDF of Mr instead of greater than Mr
CPPFLAGS += -DLF_FROM_DATA      #Use a tabluated Luminosity Function instead of a Schechter Fit
CPPFLAGS += -DRED_FRACTION      #Try to force a modeled global red fraction
endif

ifeq ($(catalog_type),"BCC")
CPPFLAGS = $(WARN) $(OPT) $(INCL) -DCATALOG $(END_FLAG)
CPPFLAGS += -DCOLORS            #Do we add colors to the galaxies?
CPPFLAGS += -DBCGS              #Add BCGs in accordance with Sarah's observations
CPPFLAGS += -DCOLORS_FROM_RELATIVE_DENSITY #Match SEDs using relative galaxy density in zbins instaed of absoloute densities
#CPPFLAGS += -DDENSPDF_IN_MAG_BIN  #Use DENSPDF of Mr instead of greater than Mr
CPPFLAGS += -DLF_FROM_DATA      #Use a tabluated Luminosity Function instead of a Schechter Fit
CPPFLAGS += -DRED_FRACTION      #Try to force a modeled global red fraction
CPPFLAGS += -DLONG64_IDS     #are particle IDs 64-bit integers?
CPPFLAGS += -DHEALPIX        #do a healpix domain-type decomposition?
CPPFLAGS += -DBCC
#OMPFLAG = -fopenmp
endif

ifeq ($(catalog_type),"SHAM_TEST")
CPPFLAGS = $(WARN) $(OPT) $(INCL) -DCATALOG $(END_FLAG)
CPPFLAGS += -DCOLORS            #Do we add colors to the galaxies?
CPPFLAGS += -DBCGS              #Add BCGs in accordance with Sarah's observations
CPPFLAGS += -DCOLORS_FROM_RELATIVE_DENSITY #Match SEDs using relative galaxy density in zbins instaed of absoloute densities
#CPPFLAGS += -DDENSPDF_IN_MAG_BIN  #Use DENSPDF of Mr instead of greater than Mr
CPPFLAGS += -DLF_FROM_DATA      #Use a tabluated Luminosity Function instead of a Schechter Fit
CPPFLAGS += -DSNAPSHOT         #Using a snapshot instead of a lightcone
CPPFLAGS += -DSHAM_TEST      #add Consuelo SHAM galaxies back to Consuelo
endif

ifeq ($(catalog_type),"snapshot")
CPPFLAGS = $(WARN) $(OPT) $(INCL) -DCATALOG $(END_FLAG)
CPPFLAGS += -DCOLORS            #Do we add colors to the galaxies?
CPPFLAGS += -DBCGS              #Add BCGs in accordance with Sarah's observations
CPPFLAGS += -DCOLORS_FROM_RELATIVE_DENSITY #Match SEDs using relative galaxy density in zbins instaed of absoloute densities
#CPPFLAGS += -DDENSPDF_IN_MAG_BIN  #Use DENSPDF of Mr instead of greater than Mr
CPPFLAGS += -DLF_FROM_DATA      #Use a tabluated Luminosity Function instead of a Schechter Fit
CPPFLAGS += -DSINGLE_Z_BIN     #we only have 1 z-bin, no need to do the first sort
CPPFLAGS += -DSNAPSHOT         #Using a snapshot instead of a lightcone
endif

ifeq ($(catalog_type),"just_colors")
CPPFLAGS = $(WARN) $(OPT) $(INCL) -DCATALOG $(END_FLAG)
CPPFLAGS += -DCOLORS            #Do we add colors to the galaxies?
CPPFLAGS += -DBCGS              #Add BCGs in accordance with Sarah's observations
CPPFLAGS += -DCOLORS_FROM_RELATIVE_DENSITY #Match SEDs using relative galaxy density in zbins instaed of absoloute densities
#CPPFLAGS += -DDENSPDF_IN_MAG_BIN  #Use DENSPDF of Mr instead of greater than Mr
###CPPFLAGS += -DLF_FROM_DATA      #Use a tabluated Luminosity Function instead of a Schechter Fit
CPPFLAGS += -DJUST_COLORS       #Use magnitudes assigned from another code
#CPPFLAGS += -DRED_FRACTION      #Try to force a modeled global red fraction
endif

ifeq ($(catalog_type),"snapshot_colors")
CPPFLAGS = $(WARN) $(OPT) $(INCL) -DCATALOG $(END_FLAG)
CPPFLAGS += -DCOLORS            #Do we add colors to the galaxies?
CPPFLAGS += -DBCGS              #Add BCGs in accordance with Sarah's observations
CPPFLAGS += -DCOLORS_FROM_RELATIVE_DENSITY #Match SEDs using relative galaxy density in zbins instaed of absoloute densities
#CPPFLAGS += -DDENSPDF_IN_MAG_BIN  #Use DENSPDF of Mr instead of greater than Mr
CPPFLAGS += -DSINGLE_Z_BIN     #we only have 1 z-bin, no need to do the first sort
CPPFLAGS += -DSNAPSHOT         #Using a snapshot instead of a lightcone

CPPFLAGS += -DJUST_COLORS       #Use magnitudes assigned from another code
endif


# Below is a list of all possible flags
#CPPFLAGS = $(WARN) $(OPT) $(INCL) -DCATALOG $(END_FLAG)
#CPPFLAGS += -DCOLORS            #Do we add colors to the galaxies?
#CPPFLAGS += -DBCGS              #Add BCGs in accordance with Sarah's observations
#CPPFLAGS += -DCOLORS_FROM_RELATIVE_DENSITY #Match SEDs using relative galaxy density in zbins instaed of absoloute densities
#CPPFLAGS += -DDENSPDF_IN_MAG_BIN  #Use DENSPDF of Mr instead of greater than Mr
#CPPFLAGS += -DSWAP             #Swap the brightest and central most halo-galaxies
#CPPFLAGS += -DMAG_LIMITED      #Make a separate file for dim galaxies that are randomly added
#CPPFLAGS += -DFULL_SKY         #Kludge for making a full sky catalog from Hubble
#CPPFLAGS += -DLF_FROM_DATA      #Use a tabluated Luminosity Function instead of a Schechter Fit
#CPPFLAGS += -DJUST_COLORS       #Use magnitudes assigned from another code
#CPPFLAGS += -DBOLSHOI  #Readin one of the bolshoi-formatted galaxy files
#CPPFLAGS += -DSINGLE_Z_BIN     #we only have 1 z-bin, no need to do the first sort
#CPPFLAGS += -DPARALLEL -i-dynamic
#CPPFLAGS += -DDEBUG_GET_GALS
#CPPFLAGS += -DDEBUG_SEDS
#CPPFLAGS += -DSNAPSHOT         #Using a snapshot instead of a lightcone
#CPPFLAGS += -DRED_FRACTION      #Try to force a modeled global red fraction
#CPPFLAGS += -DEVOLVE_DENSPDF   #readin denspdf files from multiple LasDamas snapshots
#CPPFLAGS += -DLONG64_IDS     #are particle IDs 64-bit integers?
#CPPFLAGS += -DHEALPIX        #do a healpix domain-type decomposition?
#CPPFLAGS += -DSHAM_TEST      #add Consuelo SHAM galaxies back to Consuelo
#CPPFLAGS += -DSIXPARAMS                #fit for a skewed lognormal distribution in the field objects
#CPPFLAGS += -DDEBUG
#CPPFLAGS += -DOLDPDF
#CPPFLAGS += -DDEBUG_PDF
#CPPFLAGS += -DDEBUG_PASSIGN    #Debug particle assignment
#CPPFLAGS += -DUNITTESTS
#CPPFLAGS += -DCOLORTEST
#CPPFLAGS += -DTESTREADPART
#CPPFLAGS += -DDEBUG_PIXLC
#CPPFLAGS += -DDEBUG_SHAPES
#CPPFLAGS += -DRF_TEST
#CXX = mpic++


SRCS = *.cc

EXECNAME = hv
#OBJS = hv.o particle.o  read_MGS.o read_gadget_multi.new.o color.o galaxy.o linreg.o assigngals.new.o haloocc.o neighbor_pranavkh_cleaner.o findClosest.mbusha.o read.o singleton.o mystring.o swap.o random.o  outputs.o read_galaxies.o split_volume.o cosmology.o global_vars.o ReadParameters.o ParameterDatabase.o StringDatabase.o StringDatabase.o Logger.o Utility.o simulation.o integrator.o myrand.o biniostream.o stl_util.o #read_cube.o read_warren bcgs.o
OBJS = hv.o particle.o  read_MGS.o read_gadget_multi.new.o color.o galaxy_global.o assigngals.new.o neighbor_pranavkh_cleaner.o findClosest.mbusha.o read.o singleton.o mystring.o swap.o random.o  outputs.o read_galaxies.o split_volume.o cosmology.o global_vars.o ReadParameters.o ParameterDatabase.o StringDatabase.o StringDatabase.o Logger.o Utility.o simulation.o integrator.o myrand.o biniostream.o stl_util.o owens.o kcorrect_utils.o errormodel.o shapes.o efunc.o shapeconfig.o gno.o calcparams.o healpix_utils.o galaxy.o #read_cube.o read_warren bcgs.o

EXECNAME_SPHERE = hv_sphere
OBJS_SPHERE = hv.o particle.o  read_MGS.o read_gadget_multi.new.o color.o galaxy_global.o assigngals.new.o neighbor_sphere.o findClosest.mbusha.o read.o singleton.o mystring.o swap.o random.o  outputs.o read_galaxies.o split_volume.o cosmology.o global_vars.o ReadParameters.o ParameterDatabase.o StringDatabase.o StringDatabase.o Logger.o Utility.o simulation.o integrator.o myrand.o biniostream.o stl_util.o owens.o

EXECNAME_SPLINE = hv_spline
OBJS_SPLINE = hv.o particle.o  read_MGS.o read_gadget_multi.new.o color.o galaxy_spline.o assigngals.new.o neighbor_pranavkh_cleaner.o findClosest.mbusha.o read.o singleton.o mystring.o swap.o random.o  outputs.o read_galaxies.o split_volume.o cosmology.o global_vars.o ReadParameters.o ParameterDatabase.o StringDatabase.o StringDatabase.o Logger.o Utility.o simulation.o integrator.o myrand.o biniostream.o stl_util.o owens.o denspdf.o alglibinternal.o  diffequations.o \
        linalg.o statistics.o alglibmisc.o fasttransforms.o \
        optimization.o ap.o integration.o solvers.o \
        dataanalysis.o    interpolation.o   specialfunctions.o

EXEC1 = gethalos
OBJS1 = gethalos.o particle.o  read_cube.o read_gadget_multi.o read_MGS.o getgadgethalos.o mystring.o  singleton.o read.o read_warren.o getwhids.o functions.o

EXEC2 = lens
OBJS2 = lens.o particle.o  read_cube.o functions.o mystring.o  singleton.o read.o read_warren.o kernel.o read_clusters.o

EXEC3 = pdfxi
OBJS3 = pdfxi.o particle.o read_cube.o linreg.o assigngals.new.o read.o read_warren.o functions.o singleton.o mystring.o cf_fit.o vec3.o xilike.o galaxy.o neighbor_pranavkh.o likelihood.o wp.o utils.o covariance.o chainel.o findClosest.mbusha.o read_MGS.o read_gadget_multi.o random.o

EXEC6 = get_dm
OBJS6 = get_dm.o particle.o read_cube.o mystring.o utils.o read.o read_warren.o galaxy.o functions.o singleton.o

EXEC4 = lf
OBJS4 = lf.o galaxy.o functions.o

EXEC5 = fit_wrapper
OBJS5 = fit_wrapper.o cf_fit.o linreg.o

EXEC7 = cf_wrapper
OBJS7 = cf_wrapper.o linreg.o cf.o vec3.o

EXEC8 = cov_file
OBJS8 = cov_file.o chainel.o

EXEC9 = kcorrect_utils
OBJS9 = kcorrect_utils.o cosmology.o errormodel.o

EXEC10 = shapes
OBJS10 = shapes.o efunc.o shapeconfig.o gno.o calcparams.o

EXEC11 = errormodel
OBJS11 = errormodel.o 

EXEC12 = unittests
OBJS12 = unittests.o particle.o  read_MGS.o read_gadget_multi.new.o color.o galaxy_global.o assigngals.new.o neighbor_pranavkh_cleaner.o findClosest.mbusha.o read.o singleton.o mystring.o swap.o random.o  outputs.o read_galaxies.o split_volume.o cosmology.o global_vars.o ReadParameters.o ParameterDatabase.o StringDatabase.o StringDatabase.o Logger.o Utility.o simulation.o integrator.o myrand.o biniostream.o stl_util.o owens.o kcorrect_utils.o errormodel.o shapes.o efunc.o shapeconfig.o gno.o calcparams.o healpix_utils.o galaxy.o

EXEC_TEST = test_gasdev
OBJS_TEST = test_gasdev.o


### TARGETS ###

$(EXECNAME): $(OBJS)
	$(CXX) ${OMPFLAG} -o $(EXECNAME) $^ $(CPPFLAGS)  $(LDFLAGS) $(LOADLIBES) 

$(EXECNAME_SPHERE): $(OBJS_SPHERE)
	$(CXX) -o $(EXECNAME) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXECNAME_SPLINE): $(OBJS_SPLINE)
	$(CXX) -o $(EXECNAME) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)


$(EXEC1): $(OBJS1)
	$(CXX) -o $(EXEC1) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)


$(EXEC2): $(OBJS2)
	$(CXX) -o $(EXEC2) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC3): $(OBJS3)
	$(CXX) -o $(EXEC3) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC4): $(OBJS4)
	$(CXX) -o $(EXEC4) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC5): $(OBJS5)
	$(CXX) -o $(EXEC5) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC6): $(OBJS6)
	$(CXX) -o $(EXEC6) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC7): $(OBJS7)
	$(CXX) -o $(EXEC7) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC8): $(OBJS8)
	$(CXX) -o $(EXEC8) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC9): $(OBJS9)
	$(CXX) -o $(EXEC9) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC10): $(OBJS10)
	$(CXX) -o $(EXEC10) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC11): $(OBJS11)
	$(CXX) -o $(EXEC11) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC12): $(OBJS12)
	$(CXX) -o $(EXEC12) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

$(EXEC_TEST): $(OBJS_TEST)
	$(CXX) -o $(EXEC_TEST) $^ $(CPPFLAGS) $(LDFLAGS) $(LOADLIBES)

dep:
	makedepend -- $(INCL) -- $(SRCS)

clean:
	$(RM) $(EXECNAME) $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC9) $(EXEC10) $(EXEC11) *.o core.*

cleanup:
	$(RM) $(EXECNAME) $(EXEC1) *.o *~ core.* *.bak


### DEPS ###

# DO NOT DELETE THIS LINE -- make depend depends on it

