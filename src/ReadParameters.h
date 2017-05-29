#ifndef READPARAMETERS
#define READPARAMETERS

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <strstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert> 

using namespace std;

void readParameters();

extern string simtype;
extern string simlabel;
//extern string simnum ;
extern string flabel ;
extern string datadir ;
extern string halofile ;
extern string rnn_halofile ;
extern string simulationfile ;
extern string rnnfile ;
extern string out_path ;
extern string denspdffile ;
extern string lbcgfile ;
extern string colortrdir;
extern float DECMIN ;
extern float DECMAX;
extern float RAMIN ;
extern float RAMAX;
extern float ZREDMIN;
extern float ZREDMAX ;
extern float sinDECMAX;
extern float sinDECMIN;
extern float RMIN_REAL;
extern float RMAX_REAL;
extern float REDFRACTION1;
extern float REDFRACTION2;
extern float Z_REDFRACTION1;
extern float Z_REDFRACTION2;
extern float SCATTER;
extern int EVOLVECEN;
extern float QBASELINE;
extern int REDSHIFT_FIT;
extern int GLOBAL_FIT;
extern float cm0, cm1, cm2, cm3, cm4, cmz1, cmz2, cmz3, cmz4, cm1z1, cm1z2, cm2z1, cm2z2, cm3z1, cm1z3;
extern float cs0, cs1, cs2, cs3, cs4, csz1, csz2, csz3, csz4, cs1z1, cs1z2, cs2z1, cs2z2, cs3z1, cs1z3;
extern float fm0, fm1, fm2, fm3, fm4, fmz1, fmz2, fmz3, fmz4, fm1z1, fm1z2, fm2z1, fm2z2, fm3z1, fm1z3;
extern float fs0, fs1, fs2, fs3, fs4, fsz1, fsz2, fsz3, fsz4, fs1z1, fs1z2, fs2z1, fs2z2, fs3z1, fs1z3;
extern float p0, p1, p2, p3, p4, pz1, pz2, pz3, pz4, p1z1, p1z2, p2z1, p2z2, p3z1, p1z3;
#ifdef HEALPIX
extern long nSide;
extern long PixelNum;
#endif
#ifdef BCC
extern int LCNUM;
extern string PSTR;
extern string ZSTR;
#endif
#ifdef FITS_GALAXIES
extern string ffile;
#endif
extern float minrnn, maxrnn;
#ifdef SHAM_TEST
extern string sham_file;
#endif
extern double sim_redshift;

//From Constants.h
extern float ngals ; // Was int actually
extern float PARTICLE_FRACTION;
enum ev_type{NOEV, BLAN, FABER, TIME};
extern ev_type evolution; // Was enum actually
extern float Mstar;
extern float phistar;
extern float alpha;
extern float Q;
extern float Magmin;
extern float Magmin_dens;
extern float oMagmin ;
//extern float Magmin_pdf;
extern float BCG_Mass_lim;
extern float zTol; 
extern float ColorBinSize;

//Other common variables
//extern float phi_rescale ;

//hod readin/cut related variables 
extern int read_hod;
extern float mhost_cut;
extern float rnn_cut;
extern string hodfile;

#endif
