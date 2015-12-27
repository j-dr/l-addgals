#ifndef galaxy_h
#define galaxy_h
#include <string>     //for strings
#include <iostream>   //for i/o
#include <fstream>    //for file i/o
#include <vector>     //for vectors
#include <algorithm>  //for find_if()
#include <functional> //for greater()
#include <iterator>   //for distance()
#include <cmath>      //for floor()
#include <cassert>    //for assert()

//general include files 
#include "nr.h"  //for numerical recipes
#include "integrator.h"  //in utils
#include "pi.h"
#include "simplefunc.h" //for sqr

//addgals include files 
#include "constants.h" 
#include "color.h"
#include "myrand.h"
#include "particle.h"
#include "halo.h"
#include "cosmology.h"
//#include "chainel.h"

extern Cosmology cosmo;
class Galaxy{
 public: 
  Galaxy(float mag, int id, float ldens):magnitude(mag),gid(id),d8(ldens){particle = 0; central = false; mhost = 0.; halo = 0;};
  Galaxy(float mag, int id):magnitude(mag),gid(id){particle = 0; central = false; mhost = 0.; halo = 0;};
  //redshift=0;};
  Galaxy(float dens):d8(dens){particle = 0; magnitude=0;gid=0; central = false; mhost = 0; halo = 0;};

  void Write(ofstream &file)const{
    //5D/    file<<magnitude<<" "<<particle->Ra()<<" "<<particle->Dec()<<" "<<Z()<<" "<<central<<endl;
    file<<magnitude<<" "<<particle->Ra()<<" "<<particle->Dec()<<" "<<particle->Zred()<<" "<<central<<" "<<gid<<endl;
  }
  void WriteCFinfo(ofstream &file)const{
    if(particle)
      if((particle->X()>0)&&(particle->Y()>0)&&(particle->Z()>0))
	file<<particle->X()<<" "<<particle->Y()<<" "<<particle->Z()<<" "<<magnitude<<endl;
  }

  void Print()const{cout<<magnitude<<" "<<particle->Ra()<<" "<<particle->Dec()<<" "<<d8<<endl;}

  float Dist8()const{return d8;};
  void Dist8(float dist){
    d8 = dist;
  };
  float Mhost()const{return mhost;};
  void Mhost(float mass){
    mhost = mass;
  };
  float zGal()const{return zg;};
  void zGal(float zRedGal){
    zg = zRedGal;
  };
  float Mr()const{return magnitude;};
  //  bool IsClus()const{return cluster_gal;};
  int Gid()const{return gid;};
  void P(Particle *p){
    particle = p;
  };
  Particle * P(){return particle;};
  void H(Halo *h){ halo = h;};
  Halo * H(){return halo;};
  float Z()const{return particle->Zred();};
  int Zbin()const{
#ifdef SNAPSHOT
    return (int) (floor(particle->Z()/zbsize));
#else
    return (int) (floor(Z()/zbsize));
#endif
  }
  float Ra()const{return particle->Ra();};
  float Dec()const {return particle->Dec();};
  void Mr(float newmag){magnitude=newmag;};
  float ComovDist(float angdist){
    float rr = cosmo.RofZ(Z())*angdist/180.*PI;
    //    if(normalization == SHIFT) rr = (RofZ(redshift+ZofR(roffset))-roffset)*angdist/180.*PI;
    return rr;
  };
  void Centralize(Halo *h){
    particle->position.Reset(h->X()/sim.LengthUnit(), 
		       h->Y()/sim.LengthUnit(), 
		       h->Z()/sim.LengthUnit());
    particle->velocity = h->Velocity();
    particle->SetHaloId(h->Id());
    central = true;
  }
  void DefineCentral(){
    central = true;
  }
  bool Central()const{
    return central;
  }
 private:
  bool central;
  float magnitude;
  int gid;
  float d8;
  float mhost;
  Particle * particle;
  float zg;
  Halo * halo;
};

inline float evolve_blan(float mag, float z){
  return mag+Q*(z-0.1);
  //return mag - 2.0*(1 - 1*(z-0.1))*(z-0.1);
}

inline float evolve_faber(float mag, float z){
  return mag+Q*(log10(z) + 1.0); //the 1.0 is log10(z = 0.1)
}

inline float evolve_a(float mag, float z){
  //  return mag+Q/(1+z); 
  return mag+Q*(1./(1+z) - 1./1.1); 
}

inline float evolve_mag(float mag, float z){
  if(evolution == NOEV)
    return mag;
  if(evolution == BLAN)
    return evolve_blan(mag, z);
  if(evolution == FABER)
    return evolve_faber(mag, z);
  if(evolution == TIME)
    return evolve_a(mag, z);
  return 0;
}

inline void EvolveGal(Galaxy*& pgalaxy){
  pgalaxy->Mr(evolve_mag(pgalaxy->Mr(),pgalaxy->Z()));
  /*
  if(evolution == BLAN)
    //pgalaxy->Mr(pgalaxy->Mr()+Q*(pgalaxy->Z()-0.1));
    pgalaxy->Mr(evolve_blan(pgalaxy->Mr(),pgalaxy->Z()-0.1));
  if(evolution == FABER)
    pgalaxy->Mr(evolve_faber(pgalaxy->Mr(),pgalaxy->Z()-0.1));
  */
}

inline float deevolve_blan(float mag, float z){
  return mag-Q*(z-0.1);
  //return mag + 2.0*(1 - 1*(z-0.1))*(z-0.1);
}

inline float deevolve_faber(float mag, float z){
  return mag-Q*(log10(z) + 1.0); //the 1.0 is log10(z = 0.1)
}

inline float deevolve_a(float mag, float z){
  return mag-Q*(1./(1+z) - 1./1.1); 
}

inline float deevolve_mag(float mag, float z){
  if(evolution == NOEV)
    return mag;
  if(evolution == BLAN)
    return deevolve_blan(mag, z);
  if(evolution == FABER)
    return deevolve_faber(mag, z);
  if(evolution == TIME)
    return deevolve_a(mag, z);
  return 0;
}

inline void DeEvolveGal(Galaxy*& pgalaxy){
  pgalaxy->Mr(deevolve_mag(pgalaxy->Mr(),pgalaxy->Z()));
}

//vector <Galaxy *> GetGalaxies(double vol, float, float, float, float, float, float, float);
//vector <Galaxy *> GetGalaxies(double vol, float, float, float);
//vector <Galaxy *> GetGalaxies(double vol, float phistar);
vector <Galaxy *> GetGalaxies(double vol);
//vector <Galaxy *> GetGalaxies(double vol, ChainEl chel);
//vector <Galaxy *> GetDimGalaxies(double vol, ChainEl chel);
//void GetDimGalaxies(vector <Galaxy*> &galaxies, double zmin, double zmax, double vf, ChainEl chel);
//vector <Galaxy *> GetGalaxiesNodens(double vol);
double LF(double M, double* dummy);
float ChooseMag();
void GetMags(unsigned int n, vector <double> &mags);
void GetMags(double vol, vector <double> &mags);
void ReadLFFile(void);
void ReadDMFile(void);
double LumNumberDensityInterp(double M);
double LumNumberDensity(double M);
double LumNumberDensity(double, double);
double NdensLum(double ndens);
double NdensMagnitude(double ndens);
#endif
