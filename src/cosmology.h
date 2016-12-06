#ifndef cosmology_h
#define cosmology_h

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "physicalconstants.h"

class Cosmology{
 public:
  Cosmology(){
    omega_m=0;
    set = false;
  };

  Cosmology(float omm, float s8, float h):
  omega_m(omm),omega_l(1-omm),omega_nu(0),sigma_8(s8),omega_bar(0.0),hubble(h),name(""){set=true;};
  Cosmology(float omm, float omb, float s8, float h):
  omega_m(omm),omega_l(1-omm),omega_nu(0),omega_bar(omb), sigma_8(s8),hubble(h),name(""){set=true;};
  Cosmology(float omm, float s8, float h, std::string n):
  omega_m(omm),omega_l(1-omm),omega_nu(0),sigma_8(s8),hubble(h),omega_bar(0.0),name(n){set=true;
      std::cout<<"set up cosmology"<<name<<std::endl;};
  Cosmology(float omm, float omb, float s8, float h, std::string n):
    omega_m(omm),omega_l(1-omm),omega_nu(0),omega_bar(omb), sigma_8(s8),
      hubble(h),name(n){set=true; std::cout<<"set up cosmology"<<name<<std::endl;};
    float OmegaM()const{return omega_m;};
      float OmegaL()const{return omega_l;};
      float OmegaNu()const{return omega_nu;};
      float Sigma8()const{return sigma_8;};
      float LittleH()const{return hubble;};
      float H0()const{return hubble*100.;};
      double HubbleScale()const{return cspeed/(hubble*100.);};
      std::string Name()const{return name;};
  void Print()const{
    if(set)
      std::cout<<"OmegaM:"<<omega_m
	       <<" Omega:"<<omega_l
	       <<" OmegaNu:"<<omega_nu
	       <<" OmegaB:"<<omega_bar
	       <<" Sigma8:"<<sigma_8
	       <<" HubbleConstant:"<<H0()
	       <<" name: "<<name<<" "<<std::endl;
    else
      std::cout<<"Cosmology not yet set."<<std::endl;
  }
  double RofZ(double z);
  double ZofR(double R);
  void GetZofR(float OmegaM, float OmegaL);
  void ReadZFile();
  double DistanceModulus(double z){
    //Luminosity Distance is (1+z)*R                                                
  //R is in units of Mpc/h --> change to units of 10pc.
    double rr = RofZ(z);
    //   if(normalization == SHIFT) rr = RofZ(redshift+ZofR(roffset))-roffset;     
    return 5*log10((1+z)*rr*1e5);
  };    
 private:
  bool set;
  float omega_m;
  float omega_l;
  float omega_nu;
  float omega_bar;
  float sigma_8;
  float hubble;
  std::string name;
  bool zfileread;
  std::vector <double> redshift;
  std::vector <double> losdist;
};

#endif
