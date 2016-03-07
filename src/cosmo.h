#ifndef cosmo_h
#define cosmo_h

#include <string>     //for strings
#include <iostream>   //for i/o
#include <fstream>    //for file i/o
#include <vector>     //for vectors
#include <functional> //for greater()
#include <cstdlib>    //for exit()
#include <iterator>   //for distance()
#include <cmath>      //for floor()
#include <cassert>    //for assert()
#include "constants.h" 
#include "pi.h"
#include "simplefunc.h" //for sqr
#include "integrator.h"
#include "cosmology.h"

extern Cosmology cosmo;

//Data for kcorr(z)
static std::vector <double> redshiftk;
static std::vector <double> kcorr;

//data for sig_crit(z)
static std::vector <double> redshift_sc;
static std::vector <double> sigcrit;

//data for lumnumberdensity(M)
static std::vector <double> magnitude;
static std::vector <double> lumnumdens;

//data for distancemodulus(z)
static std::vector <double> zdistmod_z;
static std::vector <double> zdistmod_dm;

void ReadKernelFile(void);
double SigCritInv(double z);

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
