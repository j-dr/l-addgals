#ifndef cosmo_h
#define cosmo_h

#include <string>     //for strings
#include <iostream>   //for i/o
#include <fstream>    //for file i/o
#include <vector>     //for vectors
//#include <algorithm>  //for find_if()
#include <functional> //for greater()
#include <cstdlib>    //for exit()
#include <iterator>   //for distance()
#include <cmath>      //for floor()
#include <cassert>    //for assert()
#include "constants.h" 
//#include "/home/risa/code/recipes/nr.h"
//#include "nr.h"        //for ran()
#include "pi.h"
#include "simplefunc.h" //for sqr
//#include "/home/risa/code/utils/integrator.h" 
#include "/afs/slac.stanford.edu/u/ki/mbusha/projects/addgals/RisaLibs/utils/integrator.h" 
#include "cosmology.h"

using namespace std;
extern Cosmology cosmo;
//data for distance(z)
//static vector <double> redshift;
//static vector <double> losdist;

//Data for kcorr(z)
static vector <double> redshiftk;
static vector <double> kcorr;

//data for sig_crit(z)
static vector <double> redshift_sc;
static vector <double> sigcrit;

//data for lumnumberdensity(M)
static vector <double> magnitude;
static vector <double> lumnumdens;

//data for distancemodulus(z)
static vector <double> zdistmod_z;
static vector <double> zdistmod_dm;

void ReadKernelFile(void);
double SigCritInv(double z);


//double ZofR(double R);
//void ReadZFile(void);
//double KofZ(double z);
//double RofZ(double z);
//void ReadKCorr(void);
/*
inline double DistanceModulus(double redshift){
  //Luminosity Distance is (1+z)*R
  //R is in units of Mpc/h --> change to units of 10pc.
   double rr = RofZ(redshift);
   //   if(normalization == SHIFT) rr = RofZ(redshift+ZofR(roffset))-roffset;
  return 5*log10((1+redshift)*rr*1e5);
};
*/

//takes kcorrected apparent magnitude
inline double AbsMag(double appmag, double redshift){return appmag-cosmo.DistanceModulus(redshift);}

enum GEOMETRY_type{OPEN, FLAT, EDS};
const GEOMETRY_type GEOMETRY = FLAT;
const double Omega_0 = 0.3;
const double Omega_Lambda_0 = 0.7;
const double Omega_nu = 0.0;
const double H0 = 100.0;
const double Hubble_scale = cspeed/H0; //Mpc  (H0 has units km/s/Mpc)

double Epeebles(double z, double* dummy);
double angular_size_distance(double z);
double transverse_comov_distance(double z);
//double angular_size_from_comov(double D, double z);
double angular_size_inarcsec_from_comov(double D, double z);
double comoving_distance_from_arcsec(double theta, double z);
double comoving_distance_from_rad(double theta, double z);


#endif
