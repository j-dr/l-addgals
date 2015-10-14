#ifndef hv2_h
#define hv2_h
#include "galaxy.h"
#include "choose.h" //galaxy recipes
#include "cosmo.h"
#include "vector3.h"
#include "color.h"
#include "singleton.h"
#include "constants.h"
#include "myrand.h"  //for rand stuff
#include "stl_util.h"  //for copy if

void cf(string datoutfn, string cfoutfn, float brightest, float dimmest);

//vector <SEDTuple> ReadCoeff();
//vector <float> GetNeighborDist(vector <Galaxy *> galaxies);
vector <float> GetNeighborDist(vector <Galaxy *> galaxies, vector <Galaxy *> points);
vector <float> GetNeighborDist1(vector <Galaxy *> galaxies);
vector <float> GetNeighborDistOrig(vector <Galaxy *> galaxies);
vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <float> &nndist, vector <GalSED> &galseds, vector <Halo *> &halos);
//vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <float> &nndist, vector <GalSED> &galseds);
vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <GalSED> &galseds);
vector <GalSED> ReadSED();
vector <GalSED> ReadDimSED();
//#endif

/*
string outpfn = out_path+"gal_pinfo.dat";
string outdfn = out_path+"gal_dinfo.dat";
string outgfn = out_path+"gal_ginfo1.dat";
string outghfn= out_path+"gal_hinfo.dat";
string outhfn = out_path+"halos.dat";
string outcfn = out_path+"gal_cfinfo.dat";
string outgzfn = out_path+"gal_zinfo.dat";
string outrfn = out_path+"gal_rnninfo.dat";
string outafn = out_path+"gal_assign.dat";
*/
string outpfn = "gal_pinfo.dat";
string outdfn = "gal_dinfo.dat";
string outgfn = "gal_ginfo1.dat";
string outghfn= "gal_hinfo.dat";
string outhfn = "halos.dat";
string outcfn = "gal_cfinfo.dat";
string outgzfn = "gal_zinfo.dat";
string outrfn = "gal_rnninfo.dat";
string outafn = "gal_assign.dat";


double ptsize =0;

void ReadColorDens();

void DeleteAndNullifyUnchosenParticle(Particle*& p){
  if(!p->IsGal()){
    delete p;
    p = 0;
  }
}


//compare particle pointers by local mass density
bool DLess(Particle * a, Particle * b)
{
  return a->Dist8() < b->Dist8(); 
}

void AbsMagFile();

//compare particle pointers by redshift
bool zLess(Particle * a, Particle * b)
{
  return a->Zred() < b->Zred(); 
}

void AbsMagFile();

#endif
