#ifndef simulation_h
#define simulation_h
#include <string>
#include <iostream>
#include "cosmology.h"
#include "box.h"
#include <cassert>
using namespace std;

//enum sim_type{ART, HV, WAR, MILLEN, GADGET2, MGS};


class Simulation{
 public:
  Simulation(void):simtype(" "), cosmology(), boxsize(0.), particlemass(0.), np(0){};
  Simulation(std::string simtype, Cosmology cosmo):simtype(simtype), cosmology(cosmo), boxsize(0.), particlemass(0.), np(0){};

  void Boxsize(float box){boxsize=box;};
  void ParticleMass(float mpart){particlemass=mpart;};
  void Np(int npart){np=npart;};
  void Label(std::string lab){label=lab;};

  float Boxsize(){
    return boxsize;
  };
  float CubeVolume(){
    return boxsize*boxsize*boxsize;
  };
  float ParticleMass(){
    return particlemass;
  };
  int Np()const{
    return np*np*np;
  };
  std::string Label()const{
    return label;
  };
  void PrintCosmology()const{
    cosmology.Print();
  }
  std::string Type()const{
    return simtype;
  }
  Cosmology SimCosmology()const{
    return cosmology;
  }
  float LengthUnit()const{
    if (simtype=="HV")
      return 3000.;
    else
      return boxsize;
  }
  float OmegaM(void){return cosmology.OmegaM();};
  float OmegaL(void){return cosmology.OmegaL();};

 private:
  std::string simtype;
  float boxsize; //in units of Mpc/h
  //float using_length; //fraction of boxsize that is used
  float particlemass; //in units of Msun/h
  //float forceresolution; //in units of kpc/h
  int np; //cuberoot of the number of particles
  std::string label;
  Cosmology cosmology;
};

#endif
