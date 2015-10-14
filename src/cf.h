#ifndef CF_H
#define CF_H

#include <string>
#include <cstdlib>  //for exit
#include <math.h>
#include <iostream.h>
#include <fstream>
#include "vec3.h"

#include <vector>
#include "simulation.h"
#include "box.h"
#include "constants.h"
#include "linreg.h"
//#include "/home/risa/code/nr2/nr.h"
//#include "/home/risa/code/recipes/nr.h"
#include "/afs/slac.stanford.edu/u/ki/mbusha/projects/addgals/RisaLibs/recipes/nr.h"

extern Simulation sim;
using namespace std;

#define QUIET

using std::cout;

enum paircounttype {ALL, NOSAMEHALO, SAMEHALO};
const paircounttype paircount = ALL;

string inpath = "";
string outpath ="";
//string infn = inpath+"points.dat";
//string outfn = outpath+"cf.dat";

//static const double cullfrac = 1.0; //keep only this fraction of the points
//static const double Lbox = sim.Boxsize()*BOXFR;
static const double a_exp = 1.000;

//long int seed = 12398745;
//float ran1(long* seed);
const double Pi = 3.14159;
const int NC=8; //number of cells used for search grid (one side of grid)

class data {
public: 
  vec3 pos;
  data* next;
public:
  data();
  int get_points_data(ifstream* infile, float brightest, float dimmest);
  int get_points_data(ifstream* infile);
  void attach(data* last[NC][NC][NC]);
};

data::data() {
  pos.x = 0.0;
  pos.y = 0.0;
  pos.z = 0.0;
  next = NULL;
}

/*
int data::get_points_data(ifstream* infile, float brightest, float dimmest) {
  double cullfrac = 1.0; //keep only this fraction of the points
  if(dimmest>=-21.5) cullfrac = 0.5;
  if(dimmest>=-21) cullfrac = 0.08;
  if(dimmest>=-20.5) cullfrac = 0.05;
  if(dimmest>=-20) cullfrac = 0.025;
  //  if(dimmest>=-20) cullfrac = 0.5;
  double x, y, z, mag;// junk;
  *infile>>x>>y>>z>>mag;//>>junk;
  float rn = ran1(&seedbox.h);
  if((rn<cullfrac)&&(mag>brightest)&&(mag<dimmest)){
    pos.x = x;
    pos.y = y;
    pos.z = z;
    if((pos.x>sim.Boxsize())||(pos.y>sim.Boxsize())||(pos.z>sim.Boxsize())
       ||(pos.x<0)||(pos.y<0)||(pos.z<0)){
      cerr<<"ERROR: out of bounds"<<pos.x<<" "<<pos.y<<" "<<pos.z<<" "<<"\tBOXSIZE="<<sim.Boxsize()<<endl;      exit(1);
    }
    return 1;
  }
  else return 0;
}
*/

int data::get_points_data(ifstream* infile) {
  double cullfrac = 1.0; //keep only this fraction of the points
  double x, y, z;
  *infile>>x>>y>>z;
  //float rn = NR::ran1(seed);
  //if(rn<cullfrac){
    pos.x = x;
    pos.y = y;
    pos.z = z;
    if((pos.x>sim.Boxsize())||(pos.y>sim.Boxsize())||(pos.z>sim.Boxsize())
       ||(pos.x<0)||(pos.y<0)||(pos.z<0)){
      cerr<<sim.Boxsize()<<" "<<BOXFR<<" "<<sim.Boxsize()<<"ERROR: out of bounds"<<pos.x<<" "<<pos.y<<" "<<pos.z<<" "<<"\tBOXSIZE="<<sim.Boxsize()<<endl;      //exit(1);
    }
    return 1;
    //}
    //else return 0;
}

//const int RBINS = 7; //number of bins for cf
//double rmin = 0.;
//double rmax = 1.;


const int RBINS = 11; //number of bins for cf
const double cellsize = sim.Boxsize()/NC; //size of PM cell

//double rmin = -0.895;
//double rmax = 1.315;

double dlogr; //bins
double max_scale;
double pairs[RBINS+1];
data* mesh[NC][NC][NC];

int load_data(const char* filename, float, float); //load data onto grid with limits, returns number of points
int load_data(const char* filename); //load data onto grid with limits, returns number of points
void compute_pairs(double* pairs);
int setup_lower_index(double x, int& il, int& ilp);
int setup_higher_index(double x, int& ih, int& ihp);
enum looptypetype {PERIODIC_L, NORMAL, PERIODIC_H};
void setup_firstloop(int pbc, int il, int ipl, 
		     int &i2, looptypetype &looptype, double &wrap);
int setup_nextloop(int pbc, int ih, int iph, 
		    int &i2, looptypetype &looptype, double &wrap);
#endif
