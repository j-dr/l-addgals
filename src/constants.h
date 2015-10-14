#ifndef constants_h
#define constants_h
#include <ctime>
//#include <string>
//#include <iostream>
//#include <fstream>
#include <math.h> //for atan
//#include "box.h"

#include "ReadParameters.h"


#ifdef LF_FROM_DATA
#include "vector3.h"
#endif

//#define NOEV_EVLN
#define BLAN_EVLN
//#define FABER_EVLN
//#define TIME_EVLN

enum dens_measure_type{TENTH, FIFTH};
const dens_measure_type dens_measure = FIFTH;//TENTH;

// this option allows you to shift the box so that
// whatever redshift it starts at is assumed to be z=0
// resulting in a lower effective value of sigma_8.
enum normalization_type{NOSHIFT, SHIFT};
const normalization_type normalization = NOSHIFT;

static const int nhpmin=12;

#ifdef LF_FROM_DATA
//The tabulated Luminosity Function -- not really constant, but doesn't change after we read it
extern vector <float> gmagbin;
extern vector <float> gdensity;
void read_lf_data(void);
#endif


// redshift interval in which nearest neighbors are counted.
static const double delta_z=0.00333564;//1000.0/cspeed;

//float zbsize = 3*delta_z;
#ifdef SNAPSHOT
static const float zbsize = 25.0;
#else
static const float zbsize = 0.02;
#endif

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
  static const int MAXSIZE=900000000;
//static const int MAXSIZE=100000000; //requested # of particle elements
//static const int MAXSIZE=400000000; //from ja
//static const int MAXSIZE=100000000; //from pow16aaa

// physical constants
static const double angle_const = 45.0/atan(1.0);
static const double d2r =  1/angle_const;//PI/180.0;


#endif
