#ifndef constants_h
#define constants_h
#include <ctime>
//#include <string>
//#include <iostream>
//#include <fstream>
#include <math.h> //for atan
//#include "box.h"

#ifdef LF_FROM_DATA
#include "vector3.h"
#endif

//#define NOEV_EVLN
//#define BLAN_EVLN
//#define FABER_EVLN
#define TIME_EVLN

enum dens_measure_type{TENTH, FIFTH};
const dens_measure_type dens_measure = FIFTH;//TENTH;

// int ngals = 31053;
// number of points to choose for the correlation function
const int ngals = 20000;

//do we want to downsample the gadget files?
//static const float PARTICLE_FRACTION = 0.0711550;  //Consuelo down to -20
//static const float PARTICLE_FRACTION = 0.182055;  //Consuelo down to -19
static const float PARTICLE_FRACTION = 1.0;  //all particles

// this option allows you to assume either
// no evolution or a Q value, set below.
// NOEV is equivalent to setting Q = 0.
enum ev_type{NOEV, BLAN, FABER, TIME};
#ifdef BLAN_EVLN
const ev_type evolution = BLAN;
#endif 
#ifdef TIME_EVLN
const ev_type evolution = TIME;
#endif
#ifdef FABER_EVLN
const ev_type evolution = FABER;
#endif
#ifdef NOEV_EVLN
const ev_type evolution = NOEV;
#endif

// this option allows you to shift the box so that
// whatever redshift it starts at is assumed to be z=0
// resulting in a lower effective value of sigma_8.
enum normalization_type{NOSHIFT, SHIFT};
const normalization_type normalization = NOSHIFT;

static const int nhpmin=12;
//luminosity function constants
#ifdef BLAN_EVLN
//static const double Mstar = -20.44;//-20.82; //BLAN evolution
//static const double Mstar = -20.73;//from Montero-Dorta
static const double Mstar = -20.41;
#endif
#ifdef NOEV_EVLN
static const double Mstar = -20.44;//-20.82; //BLAN evolution
#endif
#ifdef TIME_EVLN
static const double Mstar = -20.34; //TIME evolution 
#endif
//static const double phistar = 0.0149;//0.0147; //blanton
//static const double alpha = -1.05;//-1.00; //blanton03
//static const double phistar = 0.009*0.97366;//Montero-Dorta w/ WMAP3 volume correction
//static const double phistar = 0.0168*0.97366;
static const double phistar 
//static const double alpha = -1.23;//Montero-Dorta
static const double alpha = -1.03;//Montero-Dorta

//for matching Millennium SAM
//static const double Mstar = -21.44;//-20.82;
//static const double phistar = 0.0132;//0.0147;
//static const double alpha = -1.20;//-1.00;

#ifdef LF_FROM_DATA
//The tabulated Luminosity Function -- not really constant, but doesn't change after we read it
extern vector <float> gmagbin;
extern vector <float> gdensity;
void read_lf_data(void);
#endif

#ifdef SNAPSHOT
static const double Q = 0.; //What we want for snapshots
#else
#ifdef BLAN_EVLN
static const double Q = -1.3;//-1.62;  //from Blanton 02
#endif
#ifdef TIME_EVLN
static const double Q = 3.5;
#endif
#ifdef NOEV_EVLN
static const double Q = 0.0;
#endif
//static const double Q = 3.2;
#endif

#ifdef PARALLEL
double Magmin; 
#else
static const double Magmin = -15.9143;
#endif
static const double Magmin_dens = -19.775;
static const double oMagmin = 25.3;//25.3;
//static const double Magmin_pdf = Magmin;//652;//52;//-19.83;
static const double Magmin2 =-18;
#ifdef MAG_LIMITED
static const double oMagMin = 18.5;
#endif

static const double BCG_Mass_lim = 4e13;

//magnitude limit around which nearest 
//neighbors are counted.
//static const double miller_h = 0.70;
static const double limmag = Magmin;//-20.55-5*log10(miller_h);
// redshift interval in which nearest neighbors are counted.
static const double delta_z=0.00333564;//1000.0/cspeed;

//float zbsize = 3*delta_z;
static const float zbsize = 0.02;

static const float zcut = 0.28;
//static const float zcut = 0.129;

// simulation constants  (things that shouldn't change unless you change the simulation)
static const int BOXES=16;  //total number of boxes in the volume
static const double npbar = 0.0370370; //number density of hv particles = 1e9/LengthUnit^3
//static const double h_100 = 0.7;

// correlation function constants
static const float rmin = 0.128825; //from ja
//static const float rmin = 0.812831; //from pow16aaa
static const float rmin_bri = 0.812831;
static const float rmax = 20.4174;
static const int nbin_bri = 9;
static const int nbin = 12; //from ja
//static const int nbin = 9; //from pow16aaa

// things that shouldn't change 
// computational constants 
static int seed = 0- (time(NULL)%10000);
//static const int MAXSIZE=2000000000; //Needed for Carmen
static const int MAXSIZE=100000000; //anticipated max # of particles
//static const int MAXSIZE=400000000; //from ja
//static const int MAXSIZE=100000000; //from pow16aaa

// physical constants
static const double angle_const = 45.0/atan(1.0);
static const double d2r =  1/angle_const;//PI/180.0;

//Catalog constraings
#ifdef SNAPSHOT
static const float zTol = 100.0; 
#else
//static const float zTol = 0.01; //Accuracy for placing galaxies in z-space
static const float zTol = 0.05; //Accuracy for placing galaxies in z-space
//static const float zTol = 20.0; //Accuracy for placing galaxies in z-space
#endif
#ifdef COLORS_FROM_RELATIVE_DENSITY 
//static const float ColorBinSize = 0.05;  //old number -- reducted for parallelization
static const float ColorBinSize = 0.02;
#endif

#endif
